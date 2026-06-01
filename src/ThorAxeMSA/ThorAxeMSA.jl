module ThorAxeMSA

import CSV
import Dates
import HHsuite_jll
import HTTP
import JSON
import Scratch
import SHA
import ThorAxe

using DataFrames: DataFrame, nrow
using MIToS.MSA: AbstractMultipleSequenceAlignment, FASTA, PIRSequences, Stockholm,
                 getannotsequence, join_msas, nsequences, read_file, sequence_id,
                 sequencenames, setreference!, stringsequence, write_file
using Random: MersenneTwister
using Statistics: mean, median
using StatsBase: sample

using ..Utils: DEFAULT_PID_THRESHOLDS, ResolvedTarget, SeedSelection, ThorAxeMSAResult,
               _http_get_with_retries, _resolve_artifact_path, decode_body,
               fasta_sequence, format_pid, format_pid_dir, protein_alignment_stats,
               resolve_sequence_name, safe_rm,
               has_s_exon_annotations, s_exon_blocks_path, set_s_exon_annotations!,
               strip_ensembl_version, write_fasta, write_json, write_s_exon_blocks_tsv,
               write_text, _classify_stage_state, _pipeline_stage_dir,
               _read_stage_state, _write_stage_state

export assemble_transcript_msa,
       build_thoraxe_msa,
       compute_identity_against_reference,
       select_best_seed

const _PHYLOSOFS_RESERVED_SYMBOLS = Set([
    ' ', '\t', '\n', '\r', '\v', '\f', '\\', '*', '>', '"', '\'', ',', '-', '_',
    '/', ';', '#', '$', '.', '&', '!', '@', '[', ']'
])
const _PHYLOSOFS_FALLBACK_SYMBOLS = [c
                                     for c in vcat(collect(Char(33):Char(126)),
                                             collect(Char(0x00BC):Char(0x017F)))
                                     if !(c in _PHYLOSOFS_RESERVED_SYMBOLS)]

const _REQUIRED_ENSEMBL_FILES = (
    "sequences.fasta",
    "exonstable.tsv",
    "tree.nh",
    "ensembl_version.csv",
    "tsl.csv"
)
const _ENSEMBL_REST_BASE = "https://rest.ensembl.org"
const _ENSEMBL_JSON_HEADERS = [
    "Accept" => "application/json",
    "Content-Type" => "application/json"
]
const _BIOMART_DATASETS_URL = "https://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL"
const _BIOMART_TEXT_HEADERS = ["Accept" => "text/plain"]
const _BIOMART_DATASETS_FILE = "ENSEMBL_MART_ENSEMBL_datasets.tsv"
const _BIOMART_DATASETS_METADATA_FILE = "ENSEMBL_MART_ENSEMBL_datasets.json"
const _TRANSCRIPT_QUERY_METADATA_FILE = "iduna_transcript_query.json"

_normalize_species_name(species::Nothing) = nothing
function _normalize_species_name(species::AbstractString)
    lowercase(replace(strip(String(species)), ' ' => '_'))
end

_thoraxe_input_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_input")
_thoraxe_msa_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_msa")
function _thoraxe_candidates_dir(workdir::AbstractString)
    joinpath(_thoraxe_msa_dir(workdir), "candidates")
end
_thoraxe_logs_dir(workdir::AbstractString) = joinpath(workdir, "logs", "thoraxe")
_thoraxe_pid_runs_dir(workdir::AbstractString) = joinpath(_thoraxe_msa_dir(workdir), "runs")
function _thoraxe_input_stage_dir(workdir::AbstractString)
    _pipeline_stage_dir(workdir, "thoraxe_input")
end
function _thoraxe_msa_stage_dir(workdir::AbstractString)
    _pipeline_stage_dir(workdir, "thoraxe_msa")
end

function _thoraxe_input_outputs(input_dir::AbstractString)
    ensembl_dir = joinpath(input_dir, "Ensembl")
    outputs = Dict{String, Any}("ensembl_dir" => ensembl_dir)
    for file in _REQUIRED_ENSEMBL_FILES
        outputs[file] = joinpath(ensembl_dir, file)
    end
    outputs[_TRANSCRIPT_QUERY_METADATA_FILE] = _transcript_query_metadata_path(input_dir)
    return outputs
end

function _thoraxe_input_identity(input_dir::AbstractString, expected)
    return merge(expected, (; input_fingerprint = _bundle_fingerprint(input_dir)))
end

function _write_thoraxe_input_state(workdir::AbstractString,
        input_dir::AbstractString,
        status::Symbol,
        identity;
        action,
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing)
    return _write_stage_state(_thoraxe_input_stage_dir(workdir);
        stage = "thoraxe_input",
        stage_key = "thoraxe_input",
        status,
        identity,
        outputs = _thoraxe_input_outputs(input_dir),
        action,
        warnings,
        exception,
        workdir)
end

function _classify_thoraxe_input_stage(workdir::AbstractString,
        input_dir::AbstractString,
        identity;
        overwrite::Bool)
    overwrite && return (; reusable = false, status = :stale, warning = nothing)
    return _classify_stage_state(_thoraxe_input_stage_dir(workdir), identity,
        _thoraxe_input_outputs(input_dir); stage_label = "ThorAxe transcript_query input")
end

function _maybe_reuse_thoraxe_input(workdir::AbstractString,
        input_dir::AbstractString,
        identity,
        expected,
        cache;
        overwrite::Bool,
        manifest_message::AbstractString,
        legacy_message::AbstractString)
    overwrite && return nothing
    if cache.reusable
        @info manifest_message input_dir
        _write_thoraxe_input_state(workdir, input_dir, :done, identity; action = :reuse)
        return input_dir
    end
    if _has_valid_ensembl_bundle(input_dir) &&
       _has_matching_transcript_query_metadata(input_dir, expected)
        @info legacy_message input_dir
        _write_thoraxe_input_state(workdir, input_dir, :done,
            _thoraxe_input_identity(input_dir, expected); action = :reuse)
        return input_dir
    end
    return nothing
end

function _warn_stage_cache(cache, path::AbstractString)
    cache.warning === nothing || @warn String(cache.warning) path status=cache.status
    return nothing
end

_stage_action(cache) = cache.status === :missing ? :run : :rebuild

function _exception_summary(err)
    return (;
        type = string(typeof(err)),
        message = sprint(showerror, err))
end

function _has_valid_ensembl_bundle(bundle_root::AbstractString)::Bool
    ensembl_dir = joinpath(bundle_root, "Ensembl")
    isdir(ensembl_dir) || return false
    # ThorAxe needs all of these files before downstream steps can run.
    for file in _REQUIRED_ENSEMBL_FILES
        path = joinpath(ensembl_dir, file)
        isfile(path) && filesize(path) > 0 || return false
    end
    return true
end

function _bundle_fingerprint(bundle_root::AbstractString)
    _has_valid_ensembl_bundle(bundle_root) || return nothing
    parts = String[]
    ensembl_dir = joinpath(bundle_root, "Ensembl")
    for file in sort(collect(_REQUIRED_ENSEMBL_FILES))
        path = joinpath(ensembl_dir, file)
        push!(parts, file)
        push!(parts, string(filesize(path)))
        push!(parts, bytes2hex(SHA.sha256(read(path))))
    end
    return bytes2hex(SHA.sha256(join(parts, '\n')))
end

function _transcript_query_metadata_path(input_dir::AbstractString)
    joinpath(input_dir, _TRANSCRIPT_QUERY_METADATA_FILE)
end

function _expected_transcript_query_metadata(target::ResolvedTarget;
        specieslist::Union{Nothing, AbstractString},
        orthology::AbstractString,
        source_kind::AbstractString,
        source_path::Union{Nothing, AbstractString} = nothing,
        source_fingerprint::Union{Nothing, AbstractString} = nothing)
    return (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        species = _normalize_species_name(target.species),
        specieslist = _normalized_specieslist(specieslist),
        orthology = String(orthology),
        source_kind = String(source_kind),
        source_path = source_path === nothing ? nothing : abspath(String(source_path)),
        source_fingerprint
    )
end

function _metadata_value_matches(stored, expected)
    if expected === nothing
        return stored === nothing || stored === missing
    end
    stored === nothing && return false
    stored === missing && return false
    return String(stored) == String(expected)
end

function _metadata_matches(path::AbstractString, expected)
    isfile(path) || return false
    try
        metadata = JSON.parse(read(path, String))
        for (key, expected_value) in pairs(expected)
            string_key = String(key)
            haskey(metadata, string_key) || return false
            _metadata_value_matches(get(metadata, string_key, nothing), expected_value) ||
                return false
        end
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _write_transcript_query_metadata!(input_dir::AbstractString, expected)
    metadata = merge(expected,
        (;
            input_fingerprint = _bundle_fingerprint(input_dir),
            written_at = string(Dates.now())
        ))
    path = _transcript_query_metadata_path(input_dir)
    return write_json(path, metadata)
end

function _has_matching_transcript_query_metadata(input_dir::AbstractString, expected)
    expected_with_fingerprint = merge(expected,
        (; input_fingerprint = _bundle_fingerprint(input_dir)))
    return _metadata_matches(_transcript_query_metadata_path(input_dir),
        expected_with_fingerprint)
end

function _open_logs(f::Function, stdout_path::AbstractString, stderr_path::AbstractString)
    mkpath(dirname(stdout_path))
    mkpath(dirname(stderr_path))
    open(stdout_path, "w") do out_io
        open(stderr_path, "w") do err_io
            f(out_io, err_io)
        end
    end
end

function _run_logged_command(command::Cmd,
        stdout_log::AbstractString,
        stderr_log::AbstractString)
    _open_logs(stdout_log, stderr_log) do out_io, err_io
        process = run(pipeline(command; stdout = out_io, stderr = err_io), wait = false)
        wait(process)
        success(process) ||
            error("ThorAxe command failed: $(command). See $(stderr_log).")
    end
    return nothing
end

function _thoraxe_runner(stdout_log::AbstractString,
        stderr_log::AbstractString)
    # ThorAxe.jl builds the exact CLI command. The runner only controls logging
    # and error messages around that command.
    return command -> _run_logged_command(command, stdout_log, stderr_log)
end

function _retry_wait_seconds(attempt::Integer)
    return min(30.0, 2.0^(attempt - 1))
end

function _biomart_cache_dir()
    return Scratch.@get_scratch!("biomart_datasets")
end

function _normalized_specieslist(specieslist::Union{Nothing, AbstractString})
    specieslist === nothing && return nothing
    stripped = strip(String(specieslist))
    return isempty(stripped) ? nothing : stripped
end

function _orthology_relationships(orthology::AbstractString)
    value = strip(String(orthology))
    if value == "1:1"
        return ["ortholog_one2one"]
    elseif value == "1:n"
        return ["ortholog_one2one", "ortholog_one2many"]
    elseif value == "m:n"
        return ["ortholog_one2one", "ortholog_one2many", "ortholog_many2many"]
    end
    error("Orthology must be one of 1:1, 1:n, or m:n; got $(orthology).")
end

function _unique_nonempty_species(items)
    seen = Set{String}()
    species = String[]
    for item in items
        item === nothing && continue
        value = String(item)
        isempty(value) && continue
        value in seen && continue
        push!(seen, value)
        push!(species, value)
    end
    return species
end

function _parse_specieslist(specieslist::Union{Nothing, AbstractString})
    normalized = _normalized_specieslist(specieslist)
    normalized === nothing && return nothing

    raw_species = String[]
    # The same option accepts a comma list, a file path, or one species name.
    if occursin(',', normalized)
        append!(raw_species, split(normalized, ','))
    elseif isfile(normalized)
        append!(raw_species, readlines(normalized))
    else
        push!(raw_species, normalized)
    end

    return _unique_nonempty_species(
        _normalize_species_name(species) for species in raw_species)
end

function _specieslist_string(species::AbstractVector{<:AbstractString})
    isempty(species) && return nothing
    return join(species, ",")
end

function _fetch_biomart_datasets_text(;
        url::AbstractString = _BIOMART_DATASETS_URL,
        retries::Integer = 4,
        http_get::Function = HTTP.get)
    # BioMart can return temporary 5xx/429 responses, so fetch through retry.
    resp = _http_get_with_retries(
        url, _BIOMART_TEXT_HEADERS; retries, sleep_seconds = 1.0, http_get)
    resp.status == 200 && return decode_body(resp)
    error("BioMart datasets metadata request failed with HTTP status $(resp.status).")
end

function _parse_biomart_gene_datasets(text::AbstractString)
    datasets = Set{String}()
    for line in split(String(text), '\n')
        fields = split(strip(line), '\t')
        length(fields) >= 2 || continue
        fields[1] == "TableSet" || continue
        dataset = strip(fields[2])
        endswith(dataset, "_gene_ensembl") || continue
        push!(datasets, dataset)
    end
    return datasets
end

function _read_biomart_cache_date(metadata_path::AbstractString)
    isfile(metadata_path) || return nothing
    try
        metadata = JSON.parse(read(metadata_path, String))
        date = get(metadata, "download_date", nothing)
        date === nothing && return nothing
        return String(date)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _write_biomart_datasets_cache!(cache_dir::AbstractString,
        text::AbstractString;
        today::Dates.Date = Dates.today(),
        url::AbstractString = _BIOMART_DATASETS_URL)
    mkpath(cache_dir)
    datasets_path = joinpath(cache_dir, _BIOMART_DATASETS_FILE)
    metadata_path = joinpath(cache_dir, _BIOMART_DATASETS_METADATA_FILE)
    tmp_datasets = string(datasets_path, ".tmp")
    write_text(tmp_datasets, text)
    mv(tmp_datasets, datasets_path; force = true)
    metadata = (;
        download_date = string(today),
        download_time = string(Dates.now()),
        url,
        status = 200
    )
    tmp_metadata = string(metadata_path, ".tmp")
    write_json(tmp_metadata, metadata)
    mv(tmp_metadata, metadata_path; force = true)
    return datasets_path
end

function _read_cached_biomart_datasets(cache_dir::AbstractString)
    datasets_path = joinpath(cache_dir, _BIOMART_DATASETS_FILE)
    isfile(datasets_path) || return nothing
    datasets = _parse_biomart_gene_datasets(read(datasets_path, String))
    isempty(datasets) && return nothing
    return datasets
end

function _load_biomart_gene_datasets(;
        cache_dir::AbstractString = _biomart_cache_dir(),
        today::Dates.Date = Dates.today(),
        fetcher::Function = _fetch_biomart_datasets_text)
    datasets_path = joinpath(cache_dir, _BIOMART_DATASETS_FILE)
    metadata_path = joinpath(cache_dir, _BIOMART_DATASETS_METADATA_FILE)
    cached_date = _read_biomart_cache_date(metadata_path)
    cached_is_current = isfile(datasets_path) && cached_date == string(today)
    # Use one fresh BioMart dataset list per day to avoid repeated network hits.
    if cached_is_current
        cached = _read_cached_biomart_datasets(cache_dir)
        cached !== nothing && return (datasets = cached, warnings = String[])
    end

    try
        text = fetcher()
        datasets = _parse_biomart_gene_datasets(text)
        isempty(datasets) &&
            error("BioMart datasets metadata did not contain Ensembl Gene datasets.")
        _write_biomart_datasets_cache!(cache_dir, text; today)
        return (datasets = datasets, warnings = String[])
    catch err
        err isa InterruptException && rethrow()
        cached = _read_cached_biomart_datasets(cache_dir)
        if cached !== nothing
            stale_date = cached_date === nothing ? "unknown date" : cached_date
            warning = "BioMart datasets metadata refresh failed; using stale cache from $(stale_date). $(sprint(showerror, err))"
            return (datasets = cached, warnings = [warning])
        end
        warning = "BioMart datasets metadata refresh failed; using the unfiltered specieslist. $(sprint(showerror, err))"
        return (datasets = nothing, warnings = [warning])
    end
end

function _biomart_gene_dataset_for_species(species::AbstractString)
    normalized = _normalize_species_name(species)
    normalized === nothing && return nothing
    parts = split(normalized, '_')
    length(parts) == 2 || return nothing
    isempty(parts[1]) && return nothing
    isempty(parts[2]) && return nothing
    return string(first(parts[1]), parts[2], "_gene_ensembl")
end

function _classify_biomart_species!(
        kept::Vector{String},
        removed::Vector{String},
        unchecked::Vector{String},
        species::AbstractString,
        datasets,
        query_species::Union{Nothing, String})
    dataset = _biomart_gene_dataset_for_species(species)
    if dataset === nothing
        push!(kept, species)
        push!(unchecked, species)
        return false
    elseif dataset in datasets
        push!(kept, species)
        return false
    elseif species == query_species
        push!(kept, species)
        return true
    end
    push!(removed, species)
    return false
end

function _resolve_biomart_datasets_specieslist(target::ResolvedTarget,
        specieslist::Union{Nothing, AbstractString};
        dataset_loader::Function = _load_biomart_gene_datasets)
    requested_species = _parse_specieslist(specieslist)
    requested_species === nothing &&
        return (specieslist = _normalized_specieslist(specieslist), warnings = String[])

    loaded = dataset_loader()
    warnings = String.(loaded.warnings)
    datasets = loaded.datasets
    datasets === nothing &&
        return (specieslist = _normalized_specieslist(specieslist), warnings = warnings)

    query_species = _normalize_species_name(target.species)
    kept = String[]
    removed = String[]
    unchecked = String[]
    missing_query_dataset = false
    # Drop species that BioMart cannot serve, but keep aliases we cannot prove.
    for species in requested_species
        missing_query_dataset |= _classify_biomart_species!(
            kept, removed, unchecked, species, datasets, query_species)
    end

    if !isempty(removed)
        push!(warnings,
            "BioMart datasets filter removed species without matching Ensembl Gene datasets or recognized aliases: $(join(removed, ", ")).")
    end
    if !isempty(unchecked)
        push!(warnings,
            "BioMart datasets filter could not derive Ensembl Gene dataset names for possible species aliases: $(join(unchecked, ", ")); keeping them unchanged.")
    end
    if missing_query_dataset && query_species !== nothing
        push!(warnings,
            "BioMart datasets filter did not find an Ensembl Gene dataset for query species $(query_species); keeping it because transcript_query requires the query species.")
    end
    isempty(kept) && error("BioMart datasets filter removed all requested species.")
    return (specieslist = _specieslist_string(kept), warnings = warnings)
end

function _homology_species(data, orthology::AbstractString)
    wanted = Set(_orthology_relationships(orthology))
    species = String[]
    for item in get(data, "data", Any[])
        for homology in get(item, "homologies", Any[])
            type = get(homology, "type", nothing)
            type isa AbstractString && String(type) in wanted || continue
            target = get(homology, "target", nothing)
            target === nothing && continue
            target_species = get(target, "species", nothing)
            target_species isa AbstractString || continue
            push!(species, _normalize_species_name(target_species))
        end
    end
    return _unique_nonempty_species(species)
end

function _fetch_ensembl_homology_data(species::AbstractString,
        gene_id::AbstractString;
        retries::Integer = 4,
        sleep_seconds::Real = 1.5,
        http_get::Function = HTTP.get)
    gene_core = strip_ensembl_version(gene_id)
    url = "$(_ENSEMBL_REST_BASE)/homology/id/$(species)/$(gene_core)?type=orthologues;sequence=none"
    resp = _http_get_with_retries(
        url, _ENSEMBL_JSON_HEADERS; retries, sleep_seconds, http_get)
    resp.status == 200 && return JSON.parse(decode_body(resp))
    error("Ensembl homology specieslist filter failed for $(gene_core) in $(species) with HTTP status $(resp.status).")
end

function _fetch_ortholog_species(target::ResolvedTarget, orthology::AbstractString;
        homology_data_fetcher::Function = _fetch_ensembl_homology_data)
    species = _normalize_species_name(target.species)
    species === nothing &&
        error("Cannot run Ensembl specieslist filter because the target species is unknown.")
    data = homology_data_fetcher(species, target.ensembl_gene_id)
    return _homology_species(data, orthology)
end

function _prepend_query_species(species::AbstractVector{<:AbstractString},
        query_species::Union{Nothing, AbstractString})
    query_species === nothing && return String.(species)
    return _unique_nonempty_species(Iterators.flatten(([String(query_species)], species)))
end

function _resolve_effective_specieslist(target::ResolvedTarget,
        specieslist::Union{Nothing, AbstractString},
        orthology::AbstractString;
        homology_species_fetcher::Function = _fetch_ortholog_species)
    _orthology_relationships(orthology)
    query_species = _normalize_species_name(target.species)
    # If Ensembl homology is unavailable, fall back to the user's species list.
    ortholog_species = try
        homology_species_fetcher(target, orthology)
    catch err
        err isa InterruptException && rethrow()
        warning = "Ensembl specieslist filter failed; using the unfiltered specieslist. $(sprint(showerror, err))"
        return (specieslist = _normalized_specieslist(specieslist), warnings = [warning])
    end

    ortholog_species = _unique_nonempty_species(ortholog_species)
    if isempty(ortholog_species)
        error("Ensembl specieslist filter found no $(orthology) ortholog species for $(target.ensembl_gene_id).")
    end

    requested_species = _parse_specieslist(specieslist)
    if requested_species === nothing
        effective_species = _prepend_query_species(ortholog_species, query_species)
        return (specieslist = _specieslist_string(effective_species), warnings = String[])
    end

    allowed_targets = Set(ortholog_species)
    matching_targets = [species
                        for species in requested_species if species in allowed_targets]
    if isempty(matching_targets)
        error("Ensembl specieslist filter found no requested species with $(orthology) orthologs for $(target.ensembl_gene_id).")
    end

    removed = [species
               for species in requested_species
               if species != query_species && !(species in allowed_targets)]
    warnings = String[]
    if !isempty(removed)
        removed_text = join(removed, ", ")
        push!(warnings,
            "Removed species without $(orthology) Ensembl orthologs for $(target.ensembl_gene_id): $(removed_text).")
    end

    effective_species = _prepend_query_species(matching_targets, query_species)
    return (specieslist = _specieslist_string(effective_species), warnings = warnings)
end

function _run_transcript_query_once(gene_core::AbstractString,
        workdir::AbstractString,
        species::Union{Nothing, AbstractString},
        specieslist::Union{Nothing, AbstractString},
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        orthology::AbstractString = "1:1",
        runner::Function = _thoraxe_runner(stdout_log, stderr_log))
    _orthology_relationships(orthology)
    cd(workdir) do
        if species === nothing
            ThorAxe.transcript_query(
                gene_core; orthology = orthology, specieslist = specieslist, runner = runner)
        else
            ThorAxe.transcript_query(
                gene_core; species = species, orthology = orthology,
                specieslist = specieslist, runner = runner)
        end
    end
    return nothing
end

function _ensure_cached_thoraxe_input(source_dir::AbstractString,
        workdir::AbstractString;
        overwrite::Bool = false,
        metadata)
    dest = _thoraxe_input_dir(workdir)
    source = abspath(String(source_dir))
    _has_valid_ensembl_bundle(source) ||
        error("Cached ThorAxe input at $(source) is not a valid transcript_query bundle.")
    source_fingerprint = _bundle_fingerprint(source)
    expected = merge(metadata,
        (;
            source_path = source,
            source_fingerprint
        ))
    identity = _thoraxe_input_identity(dest, expected)
    cache = _classify_thoraxe_input_stage(workdir, dest, identity; overwrite)
    reused = _maybe_reuse_thoraxe_input(workdir, dest, identity, expected, cache;
        overwrite,
        manifest_message = "Reusing manifest-backed ThorAxe transcript_query input.",
        legacy_message = "Adopting legacy ThorAxe transcript_query input.")
    reused === nothing || return reused

    _warn_stage_cache(cache, dest)
    action = _stage_action(cache)
    try
        _write_thoraxe_input_state(workdir, dest, :running, identity; action)
        if abspath(dest) != source
            isdir(dest) && safe_rm(dest, workdir)
            mkpath(dirname(dest))
            cp(source, dest; force = true)
        end
        _has_valid_ensembl_bundle(dest) ||
            error("Copied ThorAxe input bundle at $(dest) is incomplete.")
        _write_transcript_query_metadata!(dest, expected)
        _write_thoraxe_input_state(workdir, dest, :done,
            _thoraxe_input_identity(dest, expected); action)
    catch err
        _write_thoraxe_input_state(workdir, dest, :failed, identity;
            action,
            exception = _exception_summary(err))
        rethrow()
    end
    return dest
end

function _transcript_query_attempt_action!(gene_core::AbstractString,
        query_workdir::AbstractString,
        species::Union{Nothing, AbstractString},
        active_specieslist::Union{Nothing, AbstractString},
        stdout_log::AbstractString,
        stderr_log::AbstractString,
        tmp_gene_dir::AbstractString;
        attempt::Integer,
        attempts::Integer,
        orthology::AbstractString,
        runner::Function = _thoraxe_runner(stdout_log, stderr_log))
    _run_transcript_query_once(
        gene_core, query_workdir, species, active_specieslist,
        stdout_log, stderr_log; orthology, runner)
    _has_valid_ensembl_bundle(tmp_gene_dir) && return :done
    attempt < attempts || return :failed
    @warn "transcript_query produced an invalid bundle; retrying with the same specieslist." gene=gene_core attempt
    return :retry
end

function _run_transcript_query_with_retries!(tmp_gene_dir::AbstractString,
        gene_core::AbstractString,
        query_workdir::AbstractString,
        species::Union{Nothing, AbstractString},
        initial_specieslist::Union{Nothing, AbstractString},
        attempts::Integer,
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        orthology::AbstractString,
        runner::Function = _thoraxe_runner(stdout_log, stderr_log),
        sleep_fn::Function = sleep)
    for attempt in 1:attempts
        isdir(tmp_gene_dir) && rm(tmp_gene_dir; recursive = true, force = true)
        action = _transcript_query_attempt_action!(
            gene_core, query_workdir, species, initial_specieslist,
            stdout_log, stderr_log, tmp_gene_dir;
            attempt, attempts, orthology, runner)
        action === :done && return nothing
        action === :failed && break
        sleep_fn(_retry_wait_seconds(attempt))
    end
    return nothing
end

function _run_transcript_query_stage!(target::ResolvedTarget,
        workdir::AbstractString,
        input_dir::AbstractString,
        metadata;
        specieslist::Union{Nothing, AbstractString},
        max_retries::Integer,
        orthology::AbstractString,
        transcript_query_runner::Union{Nothing, Function},
        thoraxe_runner_factory::Function,
        sleep_fn::Function,
        action)
    isdir(input_dir) && safe_rm(input_dir, workdir)
    gene_core = strip_ensembl_version(target.ensembl_gene_id)
    logs = _thoraxe_logs_dir(workdir)
    stdout_log = joinpath(logs, "transcript_query_stdout.log")
    stderr_log = joinpath(logs, "transcript_query_stderr.log")
    species = _normalize_species_name(target.species)
    attempts = max(Int(max_retries), 1)
    active_specieslist = _normalized_specieslist(specieslist)
    @info "Running ThorAxe transcript_query." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id species specieslist=active_specieslist attempts
    runner = transcript_query_runner === nothing ?
             thoraxe_runner_factory(stdout_log, stderr_log) : transcript_query_runner
    return mktempdir(workdir; prefix = "transcript_query_") do query_workdir
        tmp_gene_dir = joinpath(query_workdir, gene_core)
        _run_transcript_query_with_retries!(
            tmp_gene_dir, gene_core, query_workdir, species, active_specieslist,
            attempts, stdout_log, stderr_log;
            orthology,
            runner,
            sleep_fn)

        _has_valid_ensembl_bundle(tmp_gene_dir) ||
            error("transcript_query did not create a valid Ensembl bundle at $(tmp_gene_dir). See $(stderr_log). If failures involve the species set, try a smaller curated specieslist.")
        mv(tmp_gene_dir, input_dir; force = true)
        _write_transcript_query_metadata!(input_dir, metadata)
        _write_thoraxe_input_state(workdir, input_dir, :done,
            _thoraxe_input_identity(input_dir, metadata); action)
        return input_dir
    end
end

function _ensure_transcript_query(target::ResolvedTarget, workdir::AbstractString;
        specieslist::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        cached_input_dir::Union{Nothing, AbstractString} = nothing,
        max_retries::Integer = 2,
        orthology::AbstractString = "1:1",
        transcript_query_runner::Union{Nothing, Function} = nothing,
        thoraxe_runner_factory::Function = _thoraxe_runner,
        sleep_fn::Function = sleep)
    _orthology_relationships(orthology)
    metadata = _expected_transcript_query_metadata(target;
        specieslist,
        orthology,
        source_kind = cached_input_dir === nothing ? "transcript_query" : "cached_input")
    if cached_input_dir !== nothing
        @info "Preparing ThorAxe input from cached transcript_query bundle." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id source=cached_input_dir workdir
        return _ensure_cached_thoraxe_input(cached_input_dir, workdir; overwrite, metadata)
    end

    input_dir = _thoraxe_input_dir(workdir)
    identity = _thoraxe_input_identity(input_dir, metadata)
    cache = _classify_thoraxe_input_stage(workdir, input_dir, identity; overwrite)
    reused = _maybe_reuse_thoraxe_input(workdir, input_dir, identity, metadata, cache;
        overwrite,
        manifest_message = "Reusing manifest-backed ThorAxe transcript_query input.",
        legacy_message = "Reusing ThorAxe transcript_query input.")
    reused === nothing || return reused

    _warn_stage_cache(cache, input_dir)
    action = _stage_action(cache)
    try
        _write_thoraxe_input_state(workdir, input_dir, :running, identity; action)
        return _run_transcript_query_stage!(target, workdir, input_dir, metadata;
            specieslist,
            max_retries,
            orthology,
            transcript_query_runner,
            thoraxe_runner_factory,
            sleep_fn,
            action)
    catch err
        _write_thoraxe_input_state(workdir, input_dir, :failed, identity;
            action,
            exception = _exception_summary(err))
        rethrow()
    end
end

function _species_from_biomart_errors_file(errors_path::AbstractString)
    isfile(errors_path) || return String[]
    species = String[]
    try
        for row in CSV.File(errors_path)
            value = getproperty(row, :Species)
            value === missing && continue
            normalized = _normalize_species_name(String(value))
            normalized === nothing && continue
            push!(species, normalized)
        end
    catch err
        err isa InterruptException && rethrow()
        return String[]
    end
    return _unique_nonempty_species(species)
end

function _species_from_biomart_stderr(stderr_log::AbstractString)
    isfile(stderr_log) || return String[]
    species = String[]
    for line in eachline(stderr_log)
        found = match(r"It can not found ([A-Za-z0-9_]+) in biomart", line)
        if found !== nothing
            push!(species, _normalize_species_name(found.captures[1]))
            continue
        end
        found = match(r"Download failed for \S+ in ([A-Za-z0-9_]+)", line)
        found !== nothing && push!(species, _normalize_species_name(found.captures[1]))
    end
    return _unique_nonempty_species(species)
end

function _biomart_transcript_query_warnings(input_dir::AbstractString,
        logs_dir::AbstractString)
    errors_path = joinpath(input_dir, "Ensembl", "errors.csv")
    stderr_log = joinpath(logs_dir, "transcript_query_stderr.log")
    error_species = _species_from_biomart_errors_file(errors_path)
    stderr_species = _species_from_biomart_stderr(stderr_log)
    species = _unique_nonempty_species(Iterators.flatten((error_species, stderr_species)))
    isempty(species) && return String[]
    return [
        "BioMart transcript_query failures recorded for species: $(join(species, ", ")). See $(errors_path) and $(stderr_log)."
    ]
end

function _join_msas_consistently(
        lhs::MSAType, rhs::MSAType) where {MSAType <:
                                           AbstractMultipleSequenceAlignment}
    common = intersect(sequencenames(lhs), sequencenames(rhs))
    isempty(common) && error("Cannot join s-exon MSAs that share no sequence names.")
    return join_msas(lhs, rhs)
end

function _transcript_path_from_table(thoraxe_dir::AbstractString, transcript_id::AbstractString)
    path_table_path = joinpath(thoraxe_dir, "path_table.csv")
    isfile(path_table_path) ||
        error("Missing ThorAxe path_table.csv at $(path_table_path).")
    path_table = DataFrame(CSV.File(path_table_path))

    transcript_core = first(split(transcript_id, "."; limit = 2))
    matches = findall(eachrow(path_table)) do row
        ids = split(String(row.TranscriptIDCluster), "/")
        transcript_core in ids || String(transcript_id) in ids
    end
    isempty(matches) &&
        error("Transcript $(transcript_id) was not found in $(path_table_path).")
    length(matches) == 1 ||
        error("Transcript $(transcript_id) matched more than one ThorAxe path.")
    return String(path_table.Path[only(matches)])
end

function _transcript_exon_ids(transcript_path::AbstractString)
    exon_tokens = split(transcript_path, "/")
    return [String(exon) for exon in exon_tokens if !(exon in ("", "start", "stop"))]
end

function _transcript_exon_file(thoraxe_dir::AbstractString, exon_id::AbstractString)
    return joinpath(thoraxe_dir, "msa", "msa_s_exon_$(exon_id).fasta")
end

function _row_matches_transcript(row, transcript_id::AbstractString)
    transcript_core = strip_ensembl_version(transcript_id)
    ids = split(String(row.TranscriptIDCluster), "/")
    return transcript_core in ids || String(transcript_id) in ids
end

function _s_exon_table(thoraxe_dir::AbstractString)
    path = joinpath(thoraxe_dir, "s_exon_table.csv")
    isfile(path) || error("Missing ThorAxe s_exon_table.csv at $(path).")
    return DataFrame(CSV.File(path))
end

function _s_exon_sequence_value(value)
    ismissing(value) && return ""
    sequence = replace(String(value), "*" => "")
    return isempty(strip(sequence)) ? "" : sequence
end

function _s_exon_sequence_from_table(thoraxe_dir::AbstractString,
        exon_id::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    table = _s_exon_table(thoraxe_dir)
    gene_core = strip_ensembl_version(gene_id)
    matches = filter(eachrow(table)) do row
        String(row.S_exonID) == String(exon_id) &&
            strip_ensembl_version(String(row.GeneID)) == gene_core &&
            _row_matches_transcript(row, transcript_id)
    end
    isempty(matches) &&
        error("Could not find sequence for $(exon_id) in ThorAxe s_exon_table.csv.")
    sequences = unique(_s_exon_sequence_value(row.S_exon_Sequence) for row in matches)
    length(sequences) == 1 ||
        error("ThorAxe s_exon_table.csv has inconsistent sequences for $(exon_id).")
    return only(sequences)
end

function _read_single_sequence_msa(name::AbstractString, sequence::AbstractString)
    return mktemp() do path, io
        write(io, ">", name, "\n", sequence, "\n")
        close(io)
        read_file(path, FASTA)
    end
end

function _transcript_exon_msa(thoraxe_dir::AbstractString,
        exon_id::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    exon_file = _transcript_exon_file(thoraxe_dir, exon_id)
    isfile(exon_file) && return read_file(exon_file, FASTA)
    startswith(exon_id, "0_") ||
        error("Expected ThorAxe s-exon MSA is missing: $(exon_file).")
    sequence = _s_exon_sequence_from_table(thoraxe_dir, exon_id, gene_id, transcript_id)
    isempty(sequence) && return nothing
    return _read_single_sequence_msa(strip_ensembl_version(gene_id), sequence)
end

function _transcript_exon_segments(thoraxe_dir::AbstractString,
        exon_ids::AbstractVector{<:AbstractString},
        gene_id::AbstractString,
        transcript_id::AbstractString)
    segments = Tuple{String, AbstractMultipleSequenceAlignment}[]
    for exon_id in exon_ids
        exon_msa = _transcript_exon_msa(thoraxe_dir, exon_id, gene_id, transcript_id)
        exon_msa === nothing && continue
        push!(segments, (String(exon_id), exon_msa))
    end
    return segments
end

function _phylosofs_dir(thoraxe_dir::AbstractString)
    return joinpath(thoraxe_dir, "phylosofs")
end

function _phylosofs_transcripts_pir(thoraxe_dir::AbstractString)
    return joinpath(_phylosofs_dir(thoraxe_dir), "transcripts.pir")
end

function _phylosofs_s_exons_tsv(thoraxe_dir::AbstractString)
    return joinpath(_phylosofs_dir(thoraxe_dir), "s_exons.tsv")
end

function _has_phylosofs_outputs(thoraxe_dir::AbstractString)
    return isfile(_phylosofs_transcripts_pir(thoraxe_dir)) &&
           isfile(_phylosofs_s_exons_tsv(thoraxe_dir))
end

function _phylosofs_s_exon_code_map(thoraxe_dir::AbstractString)
    path = _phylosofs_s_exons_tsv(thoraxe_dir)
    isfile(path) || error("ThorAxe PhyloSofS s-exon map is missing: $(path).")
    code_map = Pair{Char, String}[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        fields = split(line, '\t'; limit = 2)
        length(fields) == 2 ||
            error("Invalid PhyloSofS s-exon map line in $(path): $(repr(line)).")
        symbol = only(String(fields[2]))
        push!(code_map, symbol => String(fields[1]))
    end
    isempty(code_map) && error("ThorAxe PhyloSofS s-exon map is empty: $(path).")
    return code_map
end

function _phylosofs_transcript_tokens(seq)
    tokens = String[]
    for token in split(sequence_id(seq))
        for clustered_token in split(String(token), '/')
            transcript = strip_ensembl_version(clustered_token)
            isempty(transcript) || push!(tokens, transcript)
        end
    end
    return tokens
end

function _phylosofs_transcript_sequence(thoraxe_dir::AbstractString,
        transcript_id::AbstractString)
    path = _phylosofs_transcripts_pir(thoraxe_dir)
    isfile(path) || error("ThorAxe PhyloSofS transcript PIR is missing: $(path).")
    transcript_core = strip_ensembl_version(transcript_id)
    matches = filter(read_file(path, PIRSequences)) do seq
        transcript_core in _phylosofs_transcript_tokens(seq)
    end
    isempty(matches) &&
        error("Transcript $(transcript_id) was not found in ThorAxe PhyloSofS PIR $(path).")
    length(matches) == 1 ||
        error("Transcript $(transcript_id) matched more than one entry in ThorAxe PhyloSofS PIR $(path).")
    return only(matches)
end

function _phylosofs_transcript_symbols(thoraxe_dir::AbstractString,
        transcript_id::AbstractString)
    seq = _phylosofs_transcript_sequence(thoraxe_dir, transcript_id)
    symbols = getannotsequence(seq, "Title")
    length(symbols) == length(seq) ||
        error("ThorAxe PhyloSofS annotation length for $(transcript_id) does not match its sequence length.")
    return String(symbols)
end

function _s_exon_symbol_lookup(code_map::AbstractVector{<:Pair})
    return Dict(String(s_exon_id) => symbol for (symbol, s_exon_id) in code_map)
end

function _complete_s_exon_code_map(code_map::AbstractVector{<:Pair},
        exon_ids::AbstractVector{<:AbstractString})
    complete_map = Pair{Char, String}[symbol => String(s_exon_id)
                                      for (symbol, s_exon_id) in code_map]
    mapped_exons = Set(String(s_exon_id) for (_, s_exon_id) in complete_map)
    used_symbols = Set(symbol for (symbol, _) in complete_map)
    for exon_id in exon_ids
        exon = String(exon_id)
        exon in mapped_exons && continue
        startswith(exon, "0_") ||
            error("ThorAxe PhyloSofS s-exon map has no symbol for $(exon).")
        symbol_idx = findfirst(symbol -> !(symbol in used_symbols), _PHYLOSOFS_FALLBACK_SYMBOLS)
        symbol_idx === nothing &&
            error("No unused PhyloSofS-compatible symbol is available for $(exon).")
        symbol = _PHYLOSOFS_FALLBACK_SYMBOLS[symbol_idx]
        push!(complete_map, symbol => exon)
        push!(mapped_exons, exon)
        push!(used_symbols, symbol)
    end
    return complete_map
end

function _s_exon_symbol_for_reference_residue(residue::Char, exon_symbol::Char)
    return residue == '.' ? '.' : exon_symbol
end

function _consume_phylosofs_symbol(symbol_state,
        symbols::AbstractString,
        exon_id::AbstractString,
        exon_by_symbol::AbstractDict,
        transcript_id::AbstractString)
    symbol_state === nothing &&
        error("ThorAxe PhyloSofS annotation for $(transcript_id) is shorter than the transcript MSA reference sequence.")
    observed_symbol, state = symbol_state
    pir_exon = get(exon_by_symbol, observed_symbol, nothing)
    (pir_exon === nothing || pir_exon == String(exon_id)) ||
        error("ThorAxe PhyloSofS annotation for $(transcript_id) does not match s-exon $(exon_id).")
    return observed_symbol, iterate(symbols, state)
end

function _write_projected_exon_symbols!(io::IO,
        exon_msa,
        exon_id::AbstractString,
        gene_id::AbstractString,
        symbols::AbstractString,
        exon_symbol::Char,
        exon_by_symbol::AbstractDict,
        symbol_state,
        transcript_id::AbstractString)
    reference = resolve_sequence_name(exon_msa, gene_id)
    reference === nothing &&
        error("Could not find $(gene_id) in the s-exon MSA for $(exon_id).")
    for residue in stringsequence(exon_msa, reference)
        if residue == '.' || residue == '-'
            write(io, _s_exon_symbol_for_reference_residue(residue, exon_symbol))
        else
            _,
            symbol_state = _consume_phylosofs_symbol(
                symbol_state, symbols, exon_id, exon_by_symbol, transcript_id)
            write(io, exon_symbol)
        end
    end
    return symbol_state
end

function _project_phylosofs_symbols(exon_msas::AbstractVector,
        exon_ids::AbstractVector{<:AbstractString},
        gene_id::AbstractString,
        symbols::AbstractString,
        code_map::AbstractVector{<:Pair},
        transcript_id::AbstractString)
    symbol_by_exon = _s_exon_symbol_lookup(code_map)
    exon_by_symbol = Dict(symbol => String(s_exon_id) for (symbol, s_exon_id) in code_map)
    symbol_state = iterate(symbols)
    io = IOBuffer()
    for (exon_id, exon_msa) in zip(exon_ids, exon_msas)
        exon_symbol = get(symbol_by_exon, String(exon_id), nothing)
        exon_symbol === nothing &&
            error("ThorAxe PhyloSofS s-exon map has no symbol for $(exon_id).")
        symbol_state = _write_projected_exon_symbols!(
            io, exon_msa, exon_id, gene_id, symbols, exon_symbol, exon_by_symbol,
            symbol_state, transcript_id)
    end
    symbol_state === nothing ||
        error("ThorAxe PhyloSofS annotation for $(transcript_id) is longer than the transcript MSA reference sequence.")
    return String(take!(io))
end

function _maybe_set_s_exon_annotations!(transcript_msa::AbstractMultipleSequenceAlignment,
        exon_msas::AbstractVector,
        exon_ids::AbstractVector{<:AbstractString},
        thoraxe_dir::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    _has_phylosofs_outputs(thoraxe_dir) || return transcript_msa
    code_map = _complete_s_exon_code_map(_phylosofs_s_exon_code_map(thoraxe_dir), exon_ids)
    codes = _project_phylosofs_symbols(
        exon_msas, exon_ids, gene_id,
        _phylosofs_transcript_symbols(thoraxe_dir, transcript_id),
        code_map, transcript_id)
    set_s_exon_annotations!(transcript_msa, codes, code_map)
    return transcript_msa
end

function _transcript_msa_species(thoraxe_dir::AbstractString,
        transcript_msa::AbstractMultipleSequenceAlignment)
    s_exon_path = joinpath(thoraxe_dir, "s_exon_table.csv")
    isfile(s_exon_path) || return fill("unknown", nsequences(transcript_msa))

    table = DataFrame(CSV.File(s_exon_path))
    lookup = Dict(String(row.GeneID) => String(row.Species) for row in eachrow(table))
    return [get(lookup, String(name), "unknown") for name in sequencenames(transcript_msa)]
end

function assemble_transcript_msa(thoraxe_dir::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    transcript_path = _transcript_path_from_table(thoraxe_dir, transcript_id)
    path_exon_ids = _transcript_exon_ids(transcript_path)
    isempty(path_exon_ids) && error("No s-exons were found for $(transcript_id).")

    exon_segments = _transcript_exon_segments(
        thoraxe_dir, path_exon_ids, gene_id, transcript_id)
    isempty(exon_segments) &&
        error("No non-empty s-exons were found for $(transcript_id).")
    exon_ids = first.(exon_segments)
    exon_msas = last.(exon_segments)
    transcript_msa = reduce(_join_msas_consistently, exon_msas)

    reference = resolve_sequence_name(transcript_msa, gene_id)
    reference === nothing &&
        error("Could not find $(gene_id) in the reconstructed transcript MSA.")
    setreference!(transcript_msa, reference)
    _maybe_set_s_exon_annotations!(
        transcript_msa, exon_msas, exon_ids, thoraxe_dir, gene_id, transcript_id)

    return transcript_msa, _transcript_msa_species(thoraxe_dir, transcript_msa)
end

function _read_single_fasta_sequence(path::AbstractString)
    seq = fasta_sequence(read(path, String))
    seq === nothing && error("FASTA file $(path) did not contain a sequence.")
    return seq
end

function _extract_reference_sequence(msa::AbstractMultipleSequenceAlignment,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    for id in (gene_id, transcript_id)
        name = resolve_sequence_name(msa, id)
        if name !== nothing
            sequence = replace(stringsequence(msa, name), '-' => "", '.' => "")
            return name, uppercase(String(sequence))
        end
    end
    error("Could not find a reference sequence for $(gene_id) / $(transcript_id).")
end

function _compare_protein_sequences(query_seq::AbstractString, reference_seq::AbstractString)
    protein_alignment_stats(query_seq, reference_seq)
end

function _validate_transcript_translation(target::ResolvedTarget,
        msa::AbstractMultipleSequenceAlignment;
        workdir::Union{Nothing, AbstractString} = target.workdir)
    target.uniprot_sequence_path === nothing && return String[]
    artifact_workdir = target.workdir === nothing ? workdir : target.workdir
    uniprot_sequence_path = artifact_workdir === nothing ? target.uniprot_sequence_path :
                            _resolve_artifact_path(
        target.uniprot_sequence_path, artifact_workdir)
    isfile(uniprot_sequence_path) || return String[
        "UniProt sequence file is missing; skipped ThorAxe transcript validation."]

    query_name,
    query_seq = _extract_reference_sequence(msa,
        target.ensembl_gene_id, target.transcript_id)
    reference_seq = _read_single_fasta_sequence(uniprot_sequence_path)
    stats = _compare_protein_sequences(query_seq, reference_seq)
    warnings = String[]
    if stats.insertions != 0 || stats.deletions != 0
        error(
            "ThorAxe transcript $(target.transcript_id) has indels versus UniProt $(target.uniprot_id) ",
            "(insertions=$(stats.insertions), deletions=$(stats.deletions), query=$(query_name)).")
    elseif stats.mismatches != 0
        push!(warnings,
            "ThorAxe transcript $(target.transcript_id) has $(stats.mismatches) substitutions versus UniProt $(target.uniprot_id).")
    end
    return warnings
end

function _hhsuite_query_indices(seq::AbstractString, start::Integer, stop::Integer)
    indices = Int[]
    current = Int(start)
    for residue in seq
        if residue == '-'
            push!(indices, 0)
        else
            push!(indices, current)
            current += 1
        end
    end
    current - 1 == stop || error("Could not parse HHsuite alignment positions.")
    return indices
end

function _parse_hhsuite_query_segment(line::AbstractString)
    m = match(r"\S+\s+(\d+)\s+(\S+)\s+(\d+)", line)
    m === nothing && return nothing

    start = parse(Int, m.captures[1])
    seq = m.captures[2]
    stop = parse(Int, m.captures[3])
    return (;
        cols = findfirst(seq, line), indices = _hhsuite_query_indices(seq, start, stop))
end

function _append_hhsuite_code_line!(positions, codes, line::AbstractString, cols, indices)
    occursin(r"^\s", line) || return nothing
    isempty(indices) && return nothing
    (cols === nothing || isempty(cols)) && return nothing
    append!(positions, indices)
    append!(codes, collect(line[cols]))
    return nothing
end

function _get_codes(output::AbstractString)
    in_alignment = false
    query_line = true
    cols = 1:0
    indices = Int[]
    positions = Int[]
    codes = Char[]
    for line in split(output, '\n')
        if startswith(line, "Probab=")
            in_alignment = true
            continue
        end
        in_alignment || continue
        isempty(line) && continue
        if query_line
            parsed = _parse_hhsuite_query_segment(line)
            if parsed !== nothing
                cols = parsed.cols
                indices = parsed.indices
            end
            query_line = false
        elseif startswith(line, "Confidence")
            query_line = true
        else
            _append_hhsuite_code_line!(positions, codes, line, cols, indices)
        end
    end
    return positions, codes
end

function _identity_from_codes(positions::Vector{Int}, codes::Vector{Char})
    seen = Dict{Int, Bool}()
    for (pos, code) in zip(positions, codes)
        pos == 0 && continue
        seen[pos] = get(seen, pos, false) || code == '|'
    end
    isempty(seen) && return 0.0
    return 100 * count(identity, values(seen)) / length(seen)
end

function compute_identity_against_reference(reference_fasta::AbstractString,
        sample_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        label::AbstractString = "seed")
    mktempdir() do tmp
        ref_hhm = joinpath(tmp, "reference.hhm")
        sample_hhm = joinpath(tmp, "sample.hhm")
        out_path = joinpath(tmp, "hhalign.out")
        run(pipeline(
            `$(HHsuite_jll.hhmake()) -add_cons -M 100 -i $reference_fasta -o $ref_hhm`,
            stdout = devnull, stderr = devnull))
        run(pipeline(
            `$(HHsuite_jll.hhmake()) -add_cons -M 100 -i $sample_fasta -o $sample_hhm`,
            stdout = devnull, stderr = devnull))
        run(pipeline(
            `$(HHsuite_jll.hhalign()) -glob -M 100 -i $ref_hhm -t $sample_hhm -o $out_path`,
            stdout = devnull, stderr = devnull))
        output = read(out_path, String)
        if logs_dir !== nothing
            mkpath(logs_dir)
            cp(out_path, joinpath(logs_dir, "$(label)_hhalign.out"); force = true)
        end
        positions, codes = _get_codes(output)
        return _identity_from_codes(positions, codes)
    end
end

function _candidate_pid_dir(workdir::AbstractString, pid::Real)
    joinpath(_thoraxe_candidates_dir(workdir), format_pid_dir(pid))
end

function _candidate_sample_label(sample_idx::Integer)
    sample_idx == 0 && return "full"
    return "species_subset_$(lpad(string(sample_idx), 3, '0'))"
end

function _pid_sample_paths(workdir::AbstractString, pid::Real, sample_idx::Integer)
    candidate_dir = _candidate_pid_dir(workdir, pid)
    label = _candidate_sample_label(sample_idx)
    msa_name = sample_idx == 0 ? "candidate_msa_full" : "candidate_msa_$(label)"
    sequence_name = sample_idx == 0 ? "candidate_sequences_full" :
                    "candidate_sequences_$(label)"
    species_name = sample_idx == 0 ? "candidate_species_full" :
                   "candidate_$(label)"
    return (;
        fasta_path = joinpath(candidate_dir, "$(msa_name).fasta"),
        stockholm_path = joinpath(candidate_dir, "$(msa_name).sto"),
        s_exon_blocks_tsv = joinpath(candidate_dir, "$(msa_name)_s_exon_blocks.tsv"),
        sequence_fasta = joinpath(candidate_dir, "sequences", "$(sequence_name).fasta"),
        species_file = joinpath(candidate_dir, "species", "$(species_name).txt")
    )
end

function _pid_scores_path(workdir::AbstractString, pid::Real)
    joinpath(_candidate_pid_dir(workdir, pid), "scores.csv")
end

function _pid_sample_run_root(workdir::AbstractString, pid::Real, sample_idx::Integer)
    joinpath(_thoraxe_pid_runs_dir(workdir), format_pid_dir(pid),
        _candidate_sample_label(sample_idx))
end

function _pid_sample_thoraxe_dir(workdir::AbstractString, pid::Real, sample_idx::Integer)
    joinpath(_pid_sample_run_root(workdir, pid, sample_idx), "thoraxe")
end

function _write_species_file(path::AbstractString,
        species::AbstractVector{<:AbstractString};
        overwrite::Bool = false)
    if !overwrite && isfile(path)
        return path
    end
    mkpath(dirname(path))
    open(path, "w") do io
        for item in species
            println(io, String(item))
        end
    end
    return path
end

function _write_candidate_sample_inputs(paths,
        msa::AbstractMultipleSequenceAlignment,
        species::AbstractVector{<:AbstractString},
        indices::AbstractVector{<:Integer};
        overwrite::Bool = false)
    if !overwrite && isfile(paths.sequence_fasta) && isfile(paths.species_file)
        return paths.sequence_fasta, paths.species_file
    end
    length(species) == nsequences(msa) ||
        error("Species list length does not match the MSA sequence count.")
    seq_names = collect(sequencenames(msa))
    names = String.(seq_names)
    mkpath(dirname(paths.sequence_fasta))
    write_fasta(paths.sequence_fasta,
        [(names[i], replace(stringsequence(msa, seq_names[i]), '-' => "", '.' => ""))
         for i in indices])
    _write_species_file(paths.species_file, String[species[i] for i in indices];
        overwrite = true)
    return paths.sequence_fasta, paths.species_file
end

function _write_seed_alignment_outputs(paths,
        msa::AbstractMultipleSequenceAlignment,
        pid::Real)
    write_file(paths.fasta_path, msa, FASTA)
    write_file(paths.stockholm_path, msa, Stockholm)
    write_s_exon_blocks_tsv(paths.s_exon_blocks_tsv, msa;
        alignment = "seed",
        pid = Float64(pid))
    return paths
end

function _ensure_seed_blocks_tsv(paths, pid::Real)
    isfile(paths.s_exon_blocks_tsv) && return paths.s_exon_blocks_tsv
    isfile(paths.stockholm_path) || return paths.s_exon_blocks_tsv
    msa = read_file(paths.stockholm_path, Stockholm; keepinserts = true)
    has_s_exon_annotations(msa) || return paths.s_exon_blocks_tsv
    write_s_exon_blocks_tsv(paths.s_exon_blocks_tsv, msa;
        alignment = "seed",
        pid = Float64(pid))
    return paths.s_exon_blocks_tsv
end

function _seed_has_s_exon_annotations(stockholm_path::AbstractString)
    isfile(stockholm_path) || return false
    msa = read_file(stockholm_path, Stockholm; keepinserts = true)
    return has_s_exon_annotations(msa)
end

function _reference_index(msa::AbstractMultipleSequenceAlignment,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    names = String.(sequencenames(msa))
    for id in (gene_id, transcript_id)
        name = resolve_sequence_name(msa, id)
        name === nothing && continue
        idx = findfirst(==(String(name)), names)
        idx === nothing || return idx
    end
    error("Could not find a reference sequence for $(gene_id) / $(transcript_id).")
end

function _sample_rng(seed::UInt64, sample_idx::Integer)
    mixed = xor(seed, UInt64(sample_idx) * 0xbf58476d1ce4e5b9)
    return MersenneTwister(Int(mod(mixed, UInt64(typemax(Int)))))
end

function _sample_indices(n_total::Integer, reference_idx::Integer,
        fraction::Real, rng::MersenneTwister)
    selectable = [i for i in 1:n_total if i != reference_idx]
    isempty(selectable) && return [reference_idx]
    n_keep = clamp(round(Int, Float64(fraction) * length(selectable)), 1, length(selectable))
    return vcat(reference_idx, sample(rng, selectable, n_keep; replace = false))
end

function _ensure_pid_candidate_samples(workdir::AbstractString,
        pid::Real,
        msa::AbstractMultipleSequenceAlignment,
        species::AbstractVector{<:AbstractString};
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        overwrite::Bool = false,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    reference_idx = _reference_index(msa, gene_id, transcript_id)
    @info "Preparing ThorAxe PID species samples." gene_id transcript_id pid sample_count sample_fraction
    for sample_idx in 1:sample_count
        paths = _pid_sample_paths(workdir, pid, sample_idx)
        rng = _sample_rng(sample_seed, sample_idx)
        indices = _sample_indices(nsequences(msa), reference_idx, sample_fraction, rng)
        _write_candidate_sample_inputs(paths, msa, species, indices; overwrite)
    end
    return nothing
end

function _run_kept_thoraxe_pid_msa!(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        paths,
        thoraxe_dir::AbstractString,
        path_table::AbstractString,
        pid::Real,
        specieslist::Union{Nothing, AbstractString},
        sample_idx::Integer,
        runner,
        overwrite::Bool,
        thoraxe_fn::Function = ThorAxe.thoraxe)
    if overwrite || !isfile(path_table) || !_has_phylosofs_outputs(thoraxe_dir)
        run_root = _pid_sample_run_root(workdir, pid, sample_idx)
        isdir(run_root) && safe_rm(run_root, workdir)
        mkpath(dirname(run_root))
        thoraxe_fn(
            input_dir, run_root; identity = Float64(pid), specieslist = specieslist,
            phylosofs = true, runner = runner)
    end
    pid_msa,
    _ = assemble_transcript_msa(thoraxe_dir,
        target.ensembl_gene_id, target.transcript_id)
    _write_seed_alignment_outputs(paths, pid_msa, pid)
    return paths.fasta_path, paths.stockholm_path, thoraxe_dir
end

function _run_thoraxe_pid_msa(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        specieslist::Union{Nothing, AbstractString},
        sample_idx::Integer;
        overwrite::Bool = false,
        keep_thoraxe_dir::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe)
    paths = _pid_sample_paths(workdir, pid, sample_idx)
    thoraxe_dir = _pid_sample_thoraxe_dir(workdir, pid, sample_idx)
    path_table = joinpath(thoraxe_dir, "path_table.csv")
    if !overwrite && isfile(paths.fasta_path) && isfile(paths.stockholm_path) &&
       _seed_has_s_exon_annotations(paths.stockholm_path) &&
       (!keep_thoraxe_dir || (isfile(path_table) && _has_phylosofs_outputs(thoraxe_dir)))
        @info "Reusing ThorAxe PID MSA." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid sample_idx fasta_path=paths.fasta_path stockholm_path=paths.stockholm_path
        _ensure_seed_blocks_tsv(paths, pid)
        return paths.fasta_path, paths.stockholm_path, thoraxe_dir
    end

    mkpath(_candidate_pid_dir(workdir, pid))
    pid_label = format_pid(pid)
    sample_label = "pid$(pid_label)_sample$(sample_idx)"
    stdout_log = joinpath(_thoraxe_logs_dir(workdir), "$(sample_label)_stdout.log")
    stderr_log = joinpath(_thoraxe_logs_dir(workdir), "$(sample_label)_stderr.log")
    runner = _thoraxe_runner(stdout_log, stderr_log)
    @info "Running ThorAxe PID MSA." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid sample_idx specieslist stdout_log stderr_log

    if keep_thoraxe_dir
        return _run_kept_thoraxe_pid_msa!(target, input_dir, workdir, paths,
            thoraxe_dir, path_table, pid, specieslist, sample_idx, runner, overwrite,
            thoraxe_fn)
    end

    tmp_root = joinpath(_thoraxe_msa_dir(workdir), "tmp")
    mkpath(tmp_root)
    mktempdir(tmp_root; prefix = "$(sample_label)_") do tmp
        thoraxe_fn(
            input_dir, tmp; identity = Float64(pid), specieslist = specieslist,
            phylosofs = true, runner = runner)
        pid_msa,
        _ = assemble_transcript_msa(joinpath(tmp, "thoraxe"),
            target.ensembl_gene_id, target.transcript_id)
        _write_seed_alignment_outputs(paths, pid_msa, pid)
    end
    return paths.fasta_path, paths.stockholm_path, thoraxe_dir
end

function _generate_pid_candidate(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        specieslist::Union{Nothing, AbstractString};
        overwrite::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe)
    @info "Generating ThorAxe PID candidate." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid
    paths = _pid_sample_paths(workdir, pid, 0)
    fasta_path, sto_path,
    thoraxe_dir = _run_thoraxe_pid_msa(
        target, input_dir, workdir, pid, specieslist, 0;
        overwrite, keep_thoraxe_dir = true, thoraxe_fn)
    msa,
    species = assemble_transcript_msa(thoraxe_dir,
        target.ensembl_gene_id, target.transcript_id)
    _write_seed_alignment_outputs(paths, msa, pid)
    _write_candidate_sample_inputs(
        paths, msa, species, collect(1:nsequences(msa)); overwrite)
    return (; fasta_path, stockholm_path = sto_path, thoraxe_dir,
        s_exon_blocks_tsv = paths.s_exon_blocks_tsv,
        sequence_fasta = paths.sequence_fasta, species_file = paths.species_file,
        msa, species)
end

function _candidate_msa0_validation(target::ResolvedTarget,
        msa::AbstractMultipleSequenceAlignment,
        pid::Real,
        workdir::AbstractString)
    try
        validation_warnings = _validate_transcript_translation(
            target, msa; workdir = target.workdir === nothing ? workdir : target.workdir)
        issue = isempty(validation_warnings) ? missing : join(validation_warnings, " | ")
        warnings = ["PID $(format_pid(pid)) candidate retained with warning: $(warning)"
                    for warning in validation_warnings]
        status = isempty(validation_warnings) ? "ok" : "warning"
        return (; eligible = true, status, issue, warnings)
    catch err
        err isa InterruptException && rethrow()
        issue = sprint(showerror, err)
        warning = "PID $(format_pid(pid)) excluded from seed selection: $(issue)"
        return (; eligible = false, status = "invalid_msa0", issue, warnings = [warning])
    end
end

_candidate_summary_optional(value) = value === nothing ? missing : value

function _candidate_summary_row(target::ResolvedTarget,
        candidate,
        pid::Real,
        pid_order::Integer,
        validation;
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        metadata,
        mean_identity = missing,
        median_identity = missing,
        n_samples::Integer = 0)
    return (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        pid = Float64(pid),
        pid_order,
        eligible = Bool(validation.eligible),
        selected = false,
        msa0_status = String(validation.status),
        msa0_issue = validation.issue,
        mean_identity,
        median_identity,
        n_samples = Int(n_samples),
        n_sequences_msa0 = nsequences(candidate.msa),
        pid_sample_count = Int(sample_count),
        pid_sample_fraction = Float64(sample_fraction),
        pid_sample_seed = sample_seed,
        pid_thresholds_key = metadata.pid_thresholds_key,
        effective_specieslist = _candidate_summary_optional(metadata.effective_specieslist),
        orthology = metadata.orthology,
        specieslist_filter = metadata.specieslist_filter,
        biomart_datasets_filter = metadata.biomart_datasets_filter,
        transcript_query_fingerprint = _candidate_summary_optional(
            metadata.transcript_query_fingerprint),
        selection_mode = metadata.selection_mode,
        fasta_path = candidate.fasta_path,
        stockholm_path = candidate.stockholm_path,
        sequence_fasta = candidate.sequence_fasta,
        species_file = candidate.species_file,
        scores_path = _pid_scores_path(candidate.workdir, pid)
    )
end

function _score_pid_candidate(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        pid_order::Integer,
        specieslist::Union{Nothing, AbstractString};
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        metadata,
        overwrite::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe,
        identity_fn::Function = compute_identity_against_reference)
    @info "Scoring ThorAxe PID candidate." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid pid_order sample_count
    candidate = merge(
        _generate_pid_candidate(target, input_dir, workdir, pid, specieslist;
            overwrite, thoraxe_fn),
        (; workdir))
    validation = _candidate_msa0_validation(target, candidate.msa, pid, workdir)
    if !validation.eligible
        @info "Skipping ineligible ThorAxe PID candidate." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid issue=validation.issue
        return _candidate_summary_row(
            target, candidate, pid, pid_order, validation;
            sample_count, sample_fraction, sample_seed, metadata)
    end

    if sample_count == 0
        @info "Keeping ThorAxe PID candidate without sample scoring." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid
        return _candidate_summary_row(
            target, candidate, pid, pid_order, validation;
            sample_count, sample_fraction, sample_seed, metadata)
    end

    _ensure_pid_candidate_samples(workdir, pid, candidate.msa, candidate.species;
        sample_count,
        sample_fraction,
        sample_seed,
        overwrite,
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id)

    score_rows = NamedTuple[]
    for sample_idx in 1:sample_count
        @info "Scoring ThorAxe PID sample." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid sample_idx sample_count
        sample_paths = _pid_sample_paths(workdir, pid, sample_idx)
        species_file = sample_paths.species_file
        isfile(species_file) || continue
        _run_thoraxe_pid_msa(
            target, input_dir, workdir, pid, species_file, sample_idx;
            overwrite, thoraxe_fn)
        isfile(sample_paths.fasta_path) || continue
        identity = identity_fn(
            candidate.fasta_path, sample_paths.fasta_path;
            logs_dir = joinpath(_thoraxe_logs_dir(workdir), "hhalign",
                format_pid_dir(pid)),
            label = "sample$(sample_idx)")
        push!(score_rows,
            (;
                gene_id = target.ensembl_gene_id,
                transcript_id = target.transcript_id,
                pid = Float64(pid),
                pid_order,
                sample = sample_idx,
                sample_label = _candidate_sample_label(sample_idx),
                identity,
                n_sequences_reference = nsequences(candidate.msa),
                n_sequences_sample = read_file(sample_paths.fasta_path, FASTA) |>
                                     nsequences
            ))
    end
    isempty(score_rows) && error("No PID sample scores were computed for PID $(pid).")
    CSV.write(_pid_scores_path(workdir, pid), DataFrame(score_rows))
    identities = [row.identity for row in score_rows]
    return _candidate_summary_row(
        target, candidate, pid, pid_order, validation;
        mean_identity = mean(identities),
        median_identity = median(identities),
        n_samples = length(identities),
        sample_count,
        sample_fraction,
        sample_seed,
        metadata)
end

function _summarize_candidate_scores(rows::AbstractVector{<:NamedTuple},
        path::AbstractString)
    df = DataFrame(rows)
    if :pid_order in propertynames(df)
        sort!(df, [:pid_order])
    end
    CSV.write(path, df)
    return path
end

_truthy(value) = value === true || value == 1 || lowercase(string(value)) == "true"

const _CANDIDATE_SUMMARY_STRING_TYPES = Dict(
    :pid_thresholds_key => String,
    :effective_specieslist => String,
    :orthology => String,
    :transcript_query_fingerprint => String,
    :selection_mode => String,
    :msa0_status => String,
    :msa0_issue => String,
    :stockholm_path => String,
    :fasta_path => String,
    :sequence_fasta => String,
    :species_file => String,
    :scores_path => String
)

function _candidate_summary_file(path::AbstractString)
    CSV.File(path; types = _CANDIDATE_SUMMARY_STRING_TYPES, validate = false)
end

function _candidate_summary_dataframe(path::AbstractString)
    DataFrame(_candidate_summary_file(path))
end

function _pid_thresholds_key(pid_thresholds::AbstractVector{<:Real})
    return join(repr.(Float64.(pid_thresholds)), ",")
end

function _candidate_run_metadata(input_dir::AbstractString,
        target::ResolvedTarget,
        pid_thresholds::AbstractVector{<:Real};
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        requested_sample_seed::Union{Nothing, Integer},
        effective_specieslist::Union{Nothing, AbstractString},
        orthology::AbstractString,
        specieslist_filter::Bool,
        biomart_datasets_filter::Bool)
    return (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        pid_thresholds_key = _pid_thresholds_key(pid_thresholds),
        pid_sample_count = Int(sample_count),
        pid_sample_fraction = Float64(sample_fraction),
        pid_sample_seed = sample_seed,
        requested_pid_sample_seed = requested_sample_seed,
        effective_specieslist = _normalized_specieslist(effective_specieslist),
        orthology = String(orthology),
        specieslist_filter = Bool(specieslist_filter),
        biomart_datasets_filter = Bool(biomart_datasets_filter),
        transcript_query_fingerprint = _bundle_fingerprint(input_dir),
        selection_mode = sample_count == 0 ? "all_candidates" : "sampled_selection"
    )
end

function _thoraxe_input_identity_hash(workdir::AbstractString)
    state = _read_stage_state(_thoraxe_input_stage_dir(workdir))
    state isa AbstractDict || return nothing
    value = get(state, "identity_hash", nothing)
    return value === nothing ? nothing : String(value)
end

function _canonical_thoraxe_msa_metadata(metadata)
    if metadata.requested_pid_sample_seed === nothing
        return merge(metadata, (; pid_sample_seed = nothing))
    end
    return metadata
end

function _thoraxe_msa_identity(metadata, workdir::AbstractString)
    return merge(_canonical_thoraxe_msa_metadata(metadata),
        (; thoraxe_input_identity_hash = _thoraxe_input_identity_hash(workdir)))
end

function _thoraxe_msa_outputs(summary_path::AbstractString,
        seeds::AbstractVector{SeedSelection},
        workdir::AbstractString)
    outputs = Dict{String, Any}("pid_summary" => summary_path)
    outputs["seed_stockholms"] = [_resolve_artifact_path(seed.stockholm_path, workdir)
                                  for seed in seeds]
    outputs["seed_fastas"] = [_resolve_artifact_path(seed.fasta_path, workdir)
                              for seed in seeds if seed.fasta_path !== nothing]
    return outputs
end

function _write_thoraxe_msa_state(workdir::AbstractString,
        summary_path::AbstractString,
        seeds::AbstractVector{SeedSelection},
        status::Symbol,
        identity;
        action,
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing)
    return _write_stage_state(_thoraxe_msa_stage_dir(workdir);
        stage = "thoraxe_msa",
        stage_key = "thoraxe_msa",
        status,
        identity,
        outputs = _thoraxe_msa_outputs(summary_path, seeds, workdir),
        action,
        warnings,
        exception,
        workdir)
end

function _thoraxe_msa_stage_cache(workdir::AbstractString,
        summary_path::AbstractString,
        metadata,
        stage_identity;
        overwrite::Bool)
    cache = overwrite ? (; reusable = false, status = :stale, warning = nothing) :
            _classify_stage_state(_thoraxe_msa_stage_dir(workdir), stage_identity,
        (; pid_summary = summary_path); stage_label = "ThorAxe MSA")
    has_manifest = isfile(joinpath(_thoraxe_msa_stage_dir(workdir), "stage_state.json"))
    summary_matches = !overwrite && _has_matching_candidate_summary(summary_path, metadata)
    return (; cache, has_manifest, summary_matches)
end

function _maybe_cached_thoraxe_msa(input_dir::AbstractString,
        summary_path::AbstractString,
        target::ResolvedTarget,
        workdir::AbstractString,
        metadata,
        stage_identity,
        stage_cache,
        species_filter,
        biomart_filter;
        overwrite::Bool,
        pid_sample_count::Integer,
        pid_sample_fraction::Real)
    can_try_cache = stage_cache.cache.reusable ||
                    (!overwrite && stage_cache.summary_matches)
    can_try_cache || return nothing
    cached = _cached_selected_seeds(summary_path, workdir, metadata)
    cached === nothing && return nothing
    @info "Reusing cached ThorAxe MSA candidates." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id n_seeds=length(cached.seeds) summary_path adopted_legacy=(!stage_cache.cache.reusable&&!stage_cache.has_manifest)
    result = _cached_thoraxe_msa_result(
        input_dir, summary_path, target, workdir, cached,
        species_filter, biomart_filter;
        pid_sample_count,
        pid_sample_fraction)
    _write_thoraxe_msa_state(workdir, summary_path, result.seeds, :done,
        stage_identity; action = :reuse, warnings = result.warnings)
    return result
end

function _prepare_thoraxe_msa_stage!(workdir::AbstractString,
        summary_path::AbstractString,
        stage_identity,
        stage_cache;
        overwrite::Bool)
    if stage_cache.cache.reusable
        @warn "ThorAxe MSA manifest matched, but selected seed artifacts were incomplete; rebuilding." summary_path
    elseif stage_cache.cache.warning !== nothing
        @warn String(stage_cache.cache.warning) summary_path status=stage_cache.cache.status
    end
    action = _stage_action(stage_cache.cache)
    if overwrite || (stage_cache.cache.status !== :missing && stage_cache.has_manifest)
        isdir(_thoraxe_msa_dir(workdir)) && safe_rm(_thoraxe_msa_dir(workdir), workdir)
    end
    local_artifacts_are_current = !overwrite && !stage_cache.has_manifest &&
                                  stage_cache.summary_matches
    force_pid_rerun = overwrite || !local_artifacts_are_current
    _write_thoraxe_msa_state(workdir, summary_path, SeedSelection[], :running,
        stage_identity; action)
    return (; action, local_artifacts_are_current, force_pid_rerun)
end

function _run_thoraxe_msa_stage!(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        summary_path::AbstractString,
        pid_thresholds::AbstractVector{<:Real},
        effective_specieslist::Union{Nothing, AbstractString},
        metadata,
        stage_identity,
        filters,
        prepared;
        pid_sample_count::Integer,
        pid_sample_fraction::Real,
        sample_seed::UInt64)
    score_rows = _score_pid_candidates(target, input_dir, workdir, pid_thresholds,
        effective_specieslist, metadata;
        pid_sample_count,
        pid_sample_fraction,
        sample_seed,
        overwrite = prepared.force_pid_rerun)
    _summarize_candidate_scores(score_rows, summary_path)
    @info "Selecting ThorAxe seed candidates." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id summary_path pid_sample_count
    seeds = _select_scored_candidate_seeds(summary_path, pid_sample_count)
    _mark_selected_candidates!(summary_path, seeds)
    warnings = _thoraxe_result_warnings(target, seeds, workdir, input_dir, summary_path,
        filters.species_filter, filters.biomart_filter)
    artifacts = [_seed_artifacts(seed, workdir) for seed in seeds]
    _write_thoraxe_msa_state(workdir, summary_path, seeds, :done,
        stage_identity; action = prepared.action, warnings)
    return _thoraxe_msa_result(
        input_dir, summary_path, workdir, seeds, artifacts, warnings;
        pid_sample_count,
        pid_sample_fraction,
        pid_sample_seed = sample_seed)
end

function _candidate_summary_matches(df::DataFrame, metadata)
    isempty(df) && return false
    required = (
        :gene_id,
        :transcript_id,
        :pid,
        :pid_thresholds_key,
        :pid_sample_count,
        :pid_sample_fraction,
        :pid_sample_seed,
        :effective_specieslist,
        :orthology,
        :specieslist_filter,
        :biomart_datasets_filter,
        :transcript_query_fingerprint,
        :selection_mode
    )
    all(name -> name in propertynames(df), required) || return false
    actual_pids = Float64.(df.pid)
    expected_pids = parse.(Float64, split(metadata.pid_thresholds_key, ','))
    if length(actual_pids) != length(expected_pids) ||
       Set(actual_pids) != Set(expected_pids)
        return false
    end
    return all(eachrow(df)) do row
        checks = (
            string(row.gene_id) == metadata.gene_id,
            string(row.transcript_id) == metadata.transcript_id,
            string(row.pid_thresholds_key) == metadata.pid_thresholds_key,
            Int(row.pid_sample_count) == metadata.pid_sample_count,
            Float64(row.pid_sample_fraction) == metadata.pid_sample_fraction,
            metadata.requested_pid_sample_seed === nothing ||
            UInt64(row.pid_sample_seed) == metadata.pid_sample_seed,
            _metadata_value_matches(row.effective_specieslist,
                metadata.effective_specieslist),
            string(row.orthology) == metadata.orthology,
            _truthy(row.specieslist_filter) == metadata.specieslist_filter,
            _truthy(row.biomart_datasets_filter) == metadata.biomart_datasets_filter,
            string(row.transcript_query_fingerprint) ==
            string(metadata.transcript_query_fingerprint),
            string(row.selection_mode) == metadata.selection_mode
        )
        all(checks)
    end
end

function _row_seed(row, summary_path::AbstractString)
    median_identity = ismissing(row.median_identity) ? missing :
                      Float64(row.median_identity)
    mean_identity = ismissing(row.mean_identity) ? missing : Float64(row.mean_identity)
    return SeedSelection(;
        pid = Float64(row.pid),
        median_identity,
        mean_identity,
        stockholm_path = String(row.stockholm_path),
        fasta_path = row.fasta_path === missing ? nothing : String(row.fasta_path),
        s_exon_blocks_tsv = s_exon_blocks_path(String(row.stockholm_path)),
        summary_path
    )
end

function select_best_seed(summary_path::AbstractString)
    df = _candidate_summary_dataframe(summary_path)
    isempty(df) && error("Cannot select a seed from an empty candidate summary.")
    if :eligible in propertynames(df)
        df = df[[!ismissing(value) && _truthy(value) for value in df.eligible], :]
    end
    if :median_identity in propertynames(df)
        df = df[[!ismissing(value) for value in df.median_identity], :]
    end
    if :mean_identity in propertynames(df)
        df = df[[!ismissing(value) for value in df.mean_identity], :]
    end
    isempty(df) && error("No eligible PID candidates are available for seed selection.")
    df.__row_order = 1:nrow(df)
    order_col = :pid_order in propertynames(df) ? :pid_order : :__row_order
    nseq_col = if :n_sequences_msa0 in propertynames(df)
        :n_sequences_msa0
    else
        df.__n_sequences_msa0 = zeros(Int, nrow(df))
        :__n_sequences_msa0
    end
    sort!(df, [:median_identity, :mean_identity, nseq_col, order_col];
        rev = [true, true, true, false])
    row = first(eachrow(df))
    return _row_seed(row, summary_path)
end

function _mark_selected_candidates!(summary_path::AbstractString,
        seeds::AbstractVector{SeedSelection})
    df = _candidate_summary_dataframe(summary_path)
    df.selected = falses(nrow(df))
    for seed in seeds
        selected_idx = findfirst(eachrow(df)) do row
            Float64(row.pid) == seed.pid &&
                String(row.stockholm_path) == seed.stockholm_path
        end
        selected_idx === nothing ||
            (df.selected[selected_idx] = true)
    end
    CSV.write(summary_path, df)
    return summary_path
end

function _mark_selected_candidate!(summary_path::AbstractString, seed::SeedSelection)
    _mark_selected_candidates!(summary_path, [seed])
end

function _selected_candidate_seeds(summary_path::AbstractString)
    df = _candidate_summary_dataframe(summary_path)
    :selected in propertynames(df) || return [select_best_seed(summary_path)]
    selected = df[[!ismissing(value) && _truthy(value) for value in df.selected], :]
    isempty(selected) && return SeedSelection[]
    return [_row_seed(row, summary_path) for row in eachrow(selected)]
end

function _eligible_candidate_seeds(summary_path::AbstractString)
    df = _candidate_summary_dataframe(summary_path)
    :eligible in propertynames(df) || error("Candidate summary has no eligible column.")
    eligible = df[[!ismissing(value) && _truthy(value) for value in df.eligible], :]
    isempty(eligible) && error("No eligible PID candidates are available.")
    sort!(eligible, [:pid_order])
    return [_row_seed(row, summary_path) for row in eachrow(eligible)]
end

function _normalize_pid_sample_seed(seed::Integer)::UInt64
    seed < 0 && error("pid_sample_seed must be non-negative.")
    return UInt64(seed)
end

function _validate_pid_sampling_options(sample_count::Integer, sample_fraction::Real)
    sample_count >= 0 || error("pid_sample_count must be non-negative.")
    0.0 < Float64(sample_fraction) <= 1.0 ||
        error("pid_sample_fraction must be greater than 0 and at most 1.")
    return nothing
end

function _summary_seed_value(df::DataFrame, fallback::UInt64)
    :pid_sample_seed in propertynames(df) || return fallback
    isempty(df) && return fallback
    value = df.pid_sample_seed[1]
    value === missing && return fallback
    return UInt64(value)
end

function _has_current_candidate_summary(df::DataFrame)
    required = Set([
        :gene_id,
        :transcript_id,
        :pid,
        :eligible,
        :selected,
        :msa0_status,
        :msa0_issue,
        :n_sequences_msa0,
        :pid_thresholds_key,
        :effective_specieslist,
        :orthology,
        :specieslist_filter,
        :biomart_datasets_filter,
        :transcript_query_fingerprint,
        :selection_mode,
        :sequence_fasta,
        :species_file,
        :scores_path
    ])
    names = Set(propertynames(df))
    return issubset(required, names)
end

function _candidate_summary_warnings(summary_path::AbstractString)
    isfile(summary_path) || return String[]
    df = _candidate_summary_dataframe(summary_path)
    _has_current_candidate_summary(df) || return String[]
    warnings = String[]
    for row in eachrow(df)
        issue = row.msa0_issue
        (ismissing(issue) || isempty(String(issue))) && continue
        pid_label = format_pid(Float64(row.pid))
        if !(!ismissing(row.eligible) && _truthy(row.eligible))
            push!(warnings, "PID $(pid_label) excluded from seed selection: $(issue)")
        elseif String(row.msa0_status) != "ok"
            push!(warnings, "PID $(pid_label) candidate retained with warning: $(issue)")
        end
    end
    return warnings
end

function _seed_artifacts(seed::SeedSelection, workdir::AbstractString)
    paths = _pid_sample_paths(workdir, seed.pid, 0)
    return (;
        seed_path = _resolve_artifact_path(seed.stockholm_path, workdir),
        seed_fasta = seed.fasta_path === nothing ? nothing :
                     _resolve_artifact_path(seed.fasta_path, workdir),
        s_exon_blocks_tsv = seed.s_exon_blocks_tsv === nothing ?
                            s_exon_blocks_path(
            _resolve_artifact_path(seed.stockholm_path, workdir)) :
                            _resolve_artifact_path(seed.s_exon_blocks_tsv, workdir),
        sequence_fasta = paths.sequence_fasta,
        species_file = paths.species_file,
        thoraxe_dir = _pid_sample_thoraxe_dir(workdir, seed.pid, 0)
    )
end

function _selected_artifacts_exist(seed::SeedSelection, workdir::AbstractString)
    artifacts = _seed_artifacts(seed, workdir)
    path_table = joinpath(artifacts.thoraxe_dir, "path_table.csv")
    return isfile(artifacts.seed_path) &&
           artifacts.seed_fasta !== nothing &&
           isfile(artifacts.seed_fasta) &&
           isfile(artifacts.sequence_fasta) &&
           isfile(artifacts.species_file) &&
           isfile(path_table) &&
           _seed_has_s_exon_annotations(artifacts.seed_path)
end

function _cached_selected_seeds(summary_path::AbstractString, workdir::AbstractString,
        metadata)
    isfile(summary_path) || return nothing
    df = _candidate_summary_dataframe(summary_path)
    isempty(df) && return nothing
    _has_current_candidate_summary(df) || return nothing
    _candidate_summary_matches(df, metadata) || return nothing
    seeds = _selected_candidate_seeds(summary_path)
    isempty(seeds) && return nothing
    all(seed -> _selected_artifacts_exist(seed, workdir), seeds) || return nothing
    for seed in seeds
        _ensure_seed_blocks_tsv(_pid_sample_paths(workdir, seed.pid, 0), seed.pid)
    end
    return (; seeds, sample_seed = _summary_seed_value(df, metadata.pid_sample_seed))
end

function _has_matching_candidate_summary(summary_path::AbstractString, metadata)
    isfile(summary_path) || return false
    df = _candidate_summary_dataframe(summary_path)
    isempty(df) && return false
    _has_current_candidate_summary(df) || return false
    return _candidate_summary_matches(df, metadata)
end

function _resolve_thoraxe_species_filters(target::ResolvedTarget,
        specieslist::Union{Nothing, AbstractString},
        orthology::AbstractString,
        cached_thoraxe_input_dir::Union{Nothing, AbstractString},
        specieslist_filter::Bool,
        biomart_datasets_filter::Bool;
        specieslist_resolver::Function = _resolve_effective_specieslist,
        biomart_resolver::Function = _resolve_biomart_datasets_specieslist)
    species_filter = if cached_thoraxe_input_dir === nothing && specieslist_filter
        specieslist_resolver(target, specieslist, orthology)
    else
        (specieslist = _normalized_specieslist(specieslist), warnings = String[])
    end
    biomart_filter = if cached_thoraxe_input_dir === nothing && biomart_datasets_filter
        biomart_resolver(target, species_filter.specieslist)
    else
        (specieslist = species_filter.specieslist, warnings = String[])
    end
    return (;
        species_filter,
        biomart_filter,
        effective_specieslist = biomart_filter.specieslist
    )
end

function _selected_seed_msas(seeds::AbstractVector{SeedSelection}, workdir::AbstractString)
    return [read_file(_resolve_artifact_path(seed.stockholm_path, workdir),
                Stockholm; keepinserts = true)
            for seed in seeds]
end

function _selected_translation_warnings(target::ResolvedTarget,
        selected_msas::AbstractVector,
        workdir::AbstractString)
    artifact_workdir = target.workdir === nothing ? workdir : target.workdir
    return collect(Iterators.flatten(
        _validate_transcript_translation(target, selected_msa; workdir = artifact_workdir)
    for selected_msa in selected_msas))
end

function _thoraxe_result_warnings(target::ResolvedTarget,
        seeds::AbstractVector{SeedSelection},
        workdir::AbstractString,
        input_dir::AbstractString,
        summary_path::AbstractString,
        species_filter,
        biomart_filter)
    selected_msas = _selected_seed_msas(seeds, workdir)
    return unique(vcat(species_filter.warnings,
        biomart_filter.warnings,
        _biomart_transcript_query_warnings(input_dir, _thoraxe_logs_dir(workdir)),
        _candidate_summary_warnings(summary_path),
        _selected_translation_warnings(target, selected_msas, workdir)))
end

function _thoraxe_msa_result(input_dir::AbstractString,
        summary_path::AbstractString,
        workdir::AbstractString,
        seeds::AbstractVector{SeedSelection},
        artifacts::AbstractVector,
        warnings::Vector{String};
        pid_sample_count::Integer,
        pid_sample_fraction::Real,
        pid_sample_seed)
    status = isempty(warnings) ? :ok : :warn
    return ThorAxeMSAResult(;
        input_dir,
        thoraxe_dirs = [artifact.thoraxe_dir for artifact in artifacts],
        msa_dir = _thoraxe_msa_dir(workdir),
        baseline_fastas = [String(artifact.seed_fasta) for artifact in artifacts],
        baseline_stockholms = [artifact.seed_path for artifact in artifacts],
        sequence_fastas = [artifact.sequence_fasta for artifact in artifacts],
        species_files = [artifact.species_file for artifact in artifacts],
        pid_summary = summary_path,
        seeds,
        logs_dir = _thoraxe_logs_dir(workdir),
        pid_sample_count = Int(pid_sample_count),
        pid_sample_fraction = Float64(pid_sample_fraction),
        pid_sample_seed,
        warnings,
        status
    )
end

function _cached_thoraxe_msa_result(input_dir::AbstractString,
        summary_path::AbstractString,
        target::ResolvedTarget,
        workdir::AbstractString,
        cached,
        species_filter,
        biomart_filter;
        pid_sample_count::Integer,
        pid_sample_fraction::Real)
    seeds = cached.seeds
    warnings = _thoraxe_result_warnings(
        target, seeds, workdir, input_dir, summary_path, species_filter, biomart_filter)
    artifacts = [_seed_artifacts(seed, workdir) for seed in seeds]
    return _thoraxe_msa_result(
        input_dir, summary_path, workdir, seeds, artifacts, warnings;
        pid_sample_count,
        pid_sample_fraction,
        pid_sample_seed = cached.sample_seed)
end

function _score_pid_candidates(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid_thresholds::AbstractVector{<:Real},
        effective_specieslist::Union{Nothing, AbstractString},
        metadata;
        pid_sample_count::Integer,
        pid_sample_fraction::Real,
        sample_seed::UInt64,
        overwrite::Bool,
        thoraxe_fn::Function = ThorAxe.thoraxe,
        identity_fn::Function = compute_identity_against_reference)
    score_rows = NamedTuple[]
    @info "Scoring ThorAxe PID candidates." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id n_pids=length(pid_thresholds) pid_sample_count pid_sample_fraction
    for (pid_order, pid) in enumerate(Float64.(pid_thresholds))
        push!(score_rows,
            _score_pid_candidate(
                target, input_dir, workdir, pid, pid_order, effective_specieslist;
                sample_count = Int(pid_sample_count),
                sample_fraction = Float64(pid_sample_fraction),
                sample_seed,
                metadata,
                overwrite,
                thoraxe_fn,
                identity_fn))
    end
    return score_rows
end

function _select_scored_candidate_seeds(summary_path::AbstractString,
        pid_sample_count::Integer)
    Int(pid_sample_count) == 0 && return _eligible_candidate_seeds(summary_path)
    return [select_best_seed(summary_path)]
end

function build_thoraxe_msa(target::ResolvedTarget, workdir::AbstractString;
        pid_thresholds::AbstractVector{<:Real} = DEFAULT_PID_THRESHOLDS,
        specieslist::Union{Nothing, AbstractString} = nothing,
        cached_thoraxe_input_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        orthology::AbstractString = "1:1",
        specieslist_filter::Bool = true,
        biomart_datasets_filter::Bool = true,
        transcript_query_retries::Integer = 2,
        pid_sample_count::Integer = 45,
        pid_sample_fraction::Real = 0.8,
        pid_sample_seed::Union{Nothing, Integer} = nothing)
    _orthology_relationships(orthology)
    _validate_pid_sampling_options(pid_sample_count, pid_sample_fraction)
    isempty(pid_thresholds) && error("pid_thresholds cannot be empty.")
    sample_seed = pid_sample_seed === nothing ? UInt64(rand(UInt32)) :
                  _normalize_pid_sample_seed(pid_sample_seed)
    # Species filters run only when transcript_query will create new input.
    @info "Resolving ThorAxe species filters." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id specieslist specieslist_filter biomart_datasets_filter orthology
    filters = _resolve_thoraxe_species_filters(target, specieslist, orthology,
        cached_thoraxe_input_dir, specieslist_filter, biomart_datasets_filter)
    @info "Preparing ThorAxe transcript_query input." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id effective_specieslist=filters.effective_specieslist cached_input=cached_thoraxe_input_dir
    input_dir = _ensure_transcript_query(target, workdir;
        specieslist = filters.effective_specieslist, overwrite,
        cached_input_dir = cached_thoraxe_input_dir,
        max_retries = transcript_query_retries,
        orthology)
    summary_path = joinpath(_thoraxe_msa_dir(workdir), "candidate_summary.csv")
    metadata = _candidate_run_metadata(input_dir, target, pid_thresholds;
        sample_count = Int(pid_sample_count),
        sample_fraction = Float64(pid_sample_fraction),
        sample_seed,
        requested_sample_seed = pid_sample_seed,
        effective_specieslist = filters.effective_specieslist,
        orthology,
        specieslist_filter,
        biomart_datasets_filter)
    stage_identity = _thoraxe_msa_identity(metadata, workdir)
    stage_cache = _thoraxe_msa_stage_cache(workdir, summary_path, metadata,
        stage_identity; overwrite)
    cached_result = _maybe_cached_thoraxe_msa(input_dir, summary_path, target,
        workdir, metadata, stage_identity, stage_cache, filters.species_filter,
        filters.biomart_filter;
        overwrite,
        pid_sample_count,
        pid_sample_fraction)
    cached_result === nothing || return cached_result

    prepared = _prepare_thoraxe_msa_stage!(workdir, summary_path, stage_identity,
        stage_cache; overwrite)
    try
        return _run_thoraxe_msa_stage!(target, input_dir, workdir, summary_path,
            pid_thresholds, filters.effective_specieslist, metadata, stage_identity,
            filters, prepared;
            pid_sample_count,
            pid_sample_fraction,
            sample_seed)
    catch err
        _write_thoraxe_msa_state(workdir, summary_path, SeedSelection[], :failed,
            stage_identity;
            action = prepared.action,
            exception = _exception_summary(err))
        rethrow()
    end
end

end
