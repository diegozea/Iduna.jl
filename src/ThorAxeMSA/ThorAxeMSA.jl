module ThorAxeMSA

import CSV
import Dates
import HHsuite_jll
import HTTP
import JSON3
import Scratch
import ThorAxe

using DataFrames: DataFrame
using MIToS.MSA: AbstractMultipleSequenceAlignment, FASTA, Stockholm, join_msas,
                 nsequences, read_file, sequencenames, setreference!,
                 stringsequence, write_file

using ..Utils: DEFAULT_PID_THRESHOLDS, ResolvedTarget, SeedSelection, ThorAxeMSAResult,
               _http_get_with_retries, decode_body, fasta_sequence, format_pid,
               protein_alignment_stats, resolve_sequence_name, safe_rm,
               strip_ensembl_version, write_fasta

export assemble_transcript_msa,
       build_thoraxe_msa,
       compute_identity_against_reference,
       select_best_seed

struct _CommandTimeoutError <: Exception
    command::String
    seconds::Float64
    stdout_log::String
    stderr_log::String
end

function Base.showerror(io::IO, err::_CommandTimeoutError)
    print(io,
        "ThorAxe command timed out after $(err.seconds) seconds: $(err.command). ",
        "See $(err.stdout_log) and $(err.stderr_log).")
end

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

_normalize_species_name(species::Nothing) = nothing
function _normalize_species_name(species::AbstractString)
    lowercase(replace(strip(String(species)), ' ' => '_'))
end

_thoraxe_input_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_input")
_thoraxe_output_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe")
_thoraxe_msa_dir(workdir::AbstractString) = joinpath(workdir, "thoraxe_msa")
_thoraxe_seed_dir(workdir::AbstractString) = joinpath(_thoraxe_msa_dir(workdir), "seeds")
_thoraxe_logs_dir(workdir::AbstractString) = joinpath(workdir, "logs", "thoraxe")

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

function _open_logs(f::Function, stdout_path::AbstractString, stderr_path::AbstractString)
    mkpath(dirname(stdout_path))
    mkpath(dirname(stderr_path))
    open(stdout_path, "w") do out_io
        open(stderr_path, "w") do err_io
            f(out_io, err_io)
        end
    end
end

function _normalize_timeout(timeout_seconds::Union{Nothing, Real})
    timeout_seconds === nothing && return nothing
    value = Float64(timeout_seconds)
    value > 0 || error("ThorAxe timeout values must be positive.")
    return value
end

function _wait_logged_command(process, timeout_seconds::Union{Nothing, Float64})
    if timeout_seconds === nothing
        wait(process)
        return false
    end

    deadline = time() + timeout_seconds
    # Poll so we can kill long external commands and still collect logs.
    while process_running(process)
        if time() >= deadline
            kill(process)
            wait(process)
            return true
        end
        sleep(0.25)
    end
    wait(process)
    return false
end

function _run_logged_command(command::Cmd,
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        timeout_seconds::Union{Nothing, Real} = nothing)
    timeout = _normalize_timeout(timeout_seconds)
    _open_logs(stdout_log, stderr_log) do out_io, err_io
        process = run(pipeline(command; stdout = out_io, stderr = err_io), wait = false)
        timed_out = _wait_logged_command(process, timeout)
        timed_out && throw(_CommandTimeoutError(
            string(command), timeout::Float64, String(stdout_log), String(stderr_log)))
        success(process) ||
            error("ThorAxe command failed: $(command). See $(stderr_log).")
    end
    return nothing
end

function _thoraxe_runner(stdout_log::AbstractString,
        stderr_log::AbstractString;
        timeout_seconds::Union{Nothing, Real} = nothing)
    # ThorAxe.jl builds the exact CLI command. The runner only controls logging,
    # timeouts, and error messages around that command.
    return command -> _run_logged_command(
        command, stdout_log, stderr_log; timeout_seconds = timeout_seconds)
end

function _next_timeout(timeout_seconds::Union{Nothing, Real},
        timeout_max_seconds::Union{Nothing, Real})
    timeout_seconds === nothing && return nothing
    current = _normalize_timeout(timeout_seconds)::Float64
    timeout_max_seconds === nothing && return current
    return min(2 * current, _normalize_timeout(timeout_max_seconds)::Float64)
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
        metadata = JSON3.read(read(metadata_path, String))
        date = get(metadata, :download_date, nothing)
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
    open(tmp_datasets, "w") do io
        write(io, text)
    end
    mv(tmp_datasets, datasets_path; force = true)
    metadata = (;
        download_date = string(today),
        download_time = string(Dates.now()),
        url,
        status = 200
    )
    tmp_metadata = string(metadata_path, ".tmp")
    open(tmp_metadata, "w") do io
        JSON3.pretty(io, metadata)
        println(io)
    end
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
        dataset = _biomart_gene_dataset_for_species(species)
        if dataset === nothing
            push!(kept, species)
            push!(unchecked, species)
            continue
        end
        if dataset in datasets
            push!(kept, species)
        elseif species == query_species
            push!(kept, species)
            missing_query_dataset = true
        else
            push!(removed, species)
        end
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
    resp.status == 200 && return JSON3.read(decode_body(resp))
    error("Ensembl homology specieslist filter failed for $(gene_core) in $(species) with HTTP status $(resp.status).")
end

function _fetch_ortholog_species(target::ResolvedTarget, orthology::AbstractString)
    species = _normalize_species_name(target.species)
    species === nothing &&
        error("Cannot run Ensembl specieslist filter because the target species is unknown.")
    data = _fetch_ensembl_homology_data(species, target.ensembl_gene_id)
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
        timeout_seconds::Union{Nothing, Real},
        orthology::AbstractString = "1:1",
        runner::Function = _thoraxe_runner(stdout_log, stderr_log;
            timeout_seconds = timeout_seconds))
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
        overwrite::Bool = false)
    dest = _thoraxe_input_dir(workdir)
    if !overwrite && _has_valid_ensembl_bundle(dest)
        return dest
    end

    source = abspath(String(source_dir))
    _has_valid_ensembl_bundle(source) ||
        error("Cached ThorAxe input at $(source) is not a valid transcript_query bundle.")
    if abspath(dest) != source
        isdir(dest) && safe_rm(dest, workdir)
        mkpath(dirname(dest))
        cp(source, dest; force = true)
    end
    _has_valid_ensembl_bundle(dest) ||
        error("Copied ThorAxe input bundle at $(dest) is incomplete.")
    return dest
end

function _ensure_transcript_query(target::ResolvedTarget, workdir::AbstractString;
        specieslist::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        cached_input_dir::Union{Nothing, AbstractString} = nothing,
        timeout_seconds::Union{Nothing, Real} = 180,
        timeout_max_seconds::Union{Nothing, Real} = 240,
        max_retries::Integer = 2,
        orthology::AbstractString = "1:1",
        allow_specieslist_timeout_fallback::Bool = true)
    _orthology_relationships(orthology)
    if cached_input_dir !== nothing
        return _ensure_cached_thoraxe_input(cached_input_dir, workdir; overwrite)
    end

    input_dir = _thoraxe_input_dir(workdir)
    if !overwrite && _has_valid_ensembl_bundle(input_dir)
        return input_dir
    end

    isdir(input_dir) && safe_rm(input_dir, workdir)
    gene_core = strip_ensembl_version(target.ensembl_gene_id)

    logs = _thoraxe_logs_dir(workdir)
    stdout_log = joinpath(logs, "transcript_query_stdout.log")
    stderr_log = joinpath(logs, "transcript_query_stderr.log")
    species = _normalize_species_name(target.species)

    attempts = max(Int(max_retries), 1)
    active_specieslist = _normalized_specieslist(specieslist)
    current_timeout = timeout_seconds
    # Retry transcript_query because BioMart downloads often fail transiently.
    return mktempdir(workdir; prefix = "transcript_query_") do query_workdir
        tmp_gene_dir = joinpath(query_workdir, gene_core)
        for attempt in 1:attempts
            isdir(tmp_gene_dir) && rm(tmp_gene_dir; recursive = true, force = true)
            try
                _run_transcript_query_once(
                    gene_core, query_workdir, species, active_specieslist,
                    stdout_log, stderr_log; timeout_seconds = current_timeout, orthology)
                if _has_valid_ensembl_bundle(tmp_gene_dir)
                    break
                elseif attempt < attempts
                    @warn "transcript_query produced an invalid bundle; retrying." gene=gene_core attempt
                    active_specieslist = nothing
                    sleep(_retry_wait_seconds(attempt))
                    continue
                end
            catch err
                if err isa _CommandTimeoutError && _has_valid_ensembl_bundle(tmp_gene_dir)
                    @warn "transcript_query timed out but produced a usable Ensembl bundle." gene=gene_core attempt
                    break
                end
                if err isa _CommandTimeoutError && attempt < attempts
                    if active_specieslist !== nothing && allow_specieslist_timeout_fallback
                        # A large species list can be the slow part, so retry once without it.
                        @warn "transcript_query timed out with a species list; retrying without it." gene=gene_core attempt
                        active_specieslist = nothing
                        current_timeout = timeout_max_seconds === nothing ?
                                          timeout_seconds : timeout_max_seconds
                    else
                        current_timeout = _next_timeout(current_timeout, timeout_max_seconds)
                    end
                    sleep(_retry_wait_seconds(attempt))
                    continue
                end
                rethrow(err)
            end
        end

        _has_valid_ensembl_bundle(tmp_gene_dir) ||
            error("transcript_query did not create a valid Ensembl bundle at $(tmp_gene_dir). See $(stderr_log).")
        mv(tmp_gene_dir, input_dir; force = true)
        return input_dir
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

function _ensure_baseline_thoraxe(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString;
        specieslist::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        timeout_seconds::Union{Nothing, Real} = nothing)
    thoraxe_dir = _thoraxe_output_dir(workdir)
    if !overwrite && isfile(joinpath(thoraxe_dir, "path_table.csv")) &&
       isfile(joinpath(thoraxe_dir, "s_exon_table.csv"))
        return thoraxe_dir
    end

    isdir(thoraxe_dir) && safe_rm(thoraxe_dir, workdir)
    stdout_log = joinpath(_thoraxe_logs_dir(workdir), "baseline_stdout.log")
    stderr_log = joinpath(_thoraxe_logs_dir(workdir), "baseline_stderr.log")
    runner = _thoraxe_runner(stdout_log, stderr_log; timeout_seconds = timeout_seconds)
    ThorAxe.thoraxe(input_dir, workdir; specieslist = specieslist, runner = runner)

    isfile(joinpath(thoraxe_dir, "path_table.csv")) ||
        error("ThorAxe baseline did not write path_table.csv in $(thoraxe_dir). See $(stderr_log).")
    return thoraxe_dir
end

function _join_msas_consistently(
        lhs::MSAType, rhs::MSAType) where {MSAType <:
                                           AbstractMultipleSequenceAlignment}
    common = intersect(sequencenames(lhs), sequencenames(rhs))
    isempty(common) && error("Cannot join s-exon MSAs that share no sequence names.")
    return join_msas(lhs, rhs)
end

function assemble_transcript_msa(thoraxe_dir::AbstractString,
        gene_id::AbstractString,
        transcript_id::AbstractString)
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

    transcript_path = String(path_table.Path[only(matches)])
    exon_tokens = split(replace(transcript_path, "start/" => "", "/stop" => ""), "/")
    exon_files = String[]
    # The path table lists s-exons; each one has a matching FASTA alignment.
    for exon in exon_tokens
        startswith(exon, "0_") && continue
        push!(exon_files, joinpath(thoraxe_dir, "msa", "msa_s_exon_$(exon).fasta"))
    end
    isempty(exon_files) && error("No s-exon FASTA files were found for $(transcript_id).")
    all(isfile, exon_files) || error("At least one expected s-exon MSA is missing.")

    exon_msas = [read_file(file, FASTA) for file in exon_files]
    transcript_msa = reduce(_join_msas_consistently, exon_msas)

    reference = resolve_sequence_name(transcript_msa, gene_id)
    reference === nothing &&
        error("Could not find $(gene_id) in the reconstructed transcript MSA.")
    setreference!(transcript_msa, reference)

    species = fill("unknown", nsequences(transcript_msa))
    s_exon_path = joinpath(thoraxe_dir, "s_exon_table.csv")
    if isfile(s_exon_path)
        table = DataFrame(CSV.File(s_exon_path))
        lookup = Dict(String(row.GeneID) => String(row.Species) for row in eachrow(table))
        species = [get(lookup, String(name), "unknown")
                   for name in sequencenames(transcript_msa)]
    end
    return transcript_msa, species
end

function _save_baseline_msa(workdir::AbstractString,
        msa::AbstractMultipleSequenceAlignment,
        species::AbstractVector{<:AbstractString};
        overwrite::Bool = false)
    seed_dir = _thoraxe_seed_dir(workdir)
    mkpath(seed_dir)
    fasta_path = joinpath(seed_dir, "msa_0.fasta")
    sto_path = joinpath(seed_dir, "msa_0.sto")
    sequences_path = joinpath(seed_dir, "sequences_0.fasta")
    species_path = joinpath(seed_dir, "species_0.txt")
    if !overwrite && all(isfile, (fasta_path, sto_path, sequences_path, species_path))
        return fasta_path, sto_path, sequences_path, species_path
    end

    write_file(fasta_path, msa, FASTA)
    write_file(sto_path, msa, Stockholm)
    write_fasta(sequences_path,
        [(String(name), replace(stringsequence(msa, name), '-' => ""))
         for name in sequencenames(msa)])
    open(species_path, "w") do io
        for item in species
            println(io, String(item))
        end
    end
    return fasta_path, sto_path, sequences_path, species_path
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
        msa::AbstractMultipleSequenceAlignment)
    target.uniprot_sequence_path === nothing && return String[]
    isfile(target.uniprot_sequence_path) || return String[
        "UniProt sequence file is missing; skipped ThorAxe transcript validation."]

    query_name,
    query_seq = _extract_reference_sequence(msa,
        target.ensembl_gene_id, target.transcript_id)
    reference_seq = _read_single_fasta_sequence(target.uniprot_sequence_path)
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
            m = match(r"\S+\s+(\d+)\s+(\S+)\s+(\d+)", line)
            if m !== nothing
                start = parse(Int, m.captures[1])
                seq = m.captures[2]
                stop = parse(Int, m.captures[3])
                cols = findfirst(seq, line)
                indices = Int[]
                for residue in seq
                    if residue == '-'
                        push!(indices, 0)
                    else
                        push!(indices, start)
                        start += 1
                    end
                end
                start - 1 == stop || error("Could not parse HHsuite alignment positions.")
            end
            query_line = false
        elseif startswith(line, "Confidence")
            query_line = true
        elseif occursin(r"^\s", line) && !isempty(indices) && !isempty(cols)
            append!(positions, indices)
            append!(codes, collect(line[cols]))
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

function _generate_pid_seed(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        species_file::AbstractString;
        overwrite::Bool = false,
        timeout_seconds::Union{Nothing, Real} = nothing)
    seed_dir = _thoraxe_seed_dir(workdir)
    pid_label = format_pid(pid)
    fasta_path = joinpath(seed_dir, "thoraxe_pid$(pid_label)_msa_0.fasta")
    sto_path = joinpath(seed_dir, "thoraxe_pid$(pid_label)_msa_0.sto")
    if !overwrite && isfile(fasta_path) && isfile(sto_path)
        return fasta_path, sto_path
    end

    tmp_root = joinpath(_thoraxe_msa_dir(workdir), "tmp")
    mkpath(tmp_root)
    mktempdir(tmp_root; prefix = "thoraxe_pid$(pid_label)_") do tmp
        stdout_log = joinpath(_thoraxe_logs_dir(workdir), "pid$(pid_label)_stdout.log")
        stderr_log = joinpath(_thoraxe_logs_dir(workdir), "pid$(pid_label)_stderr.log")
        runner = _thoraxe_runner(stdout_log, stderr_log; timeout_seconds = timeout_seconds)
        ThorAxe.thoraxe(
            input_dir, tmp; identity = Float64(pid), specieslist = species_file,
            runner = runner)
        pid_msa,
        _ = assemble_transcript_msa(joinpath(tmp, "thoraxe"),
            target.ensembl_gene_id, target.transcript_id)
        write_file(fasta_path, pid_msa, FASTA)
        write_file(sto_path, pid_msa, Stockholm)
    end
    return fasta_path, sto_path
end

function _summarize_pid_scores(rows::Vector{NamedTuple}, path::AbstractString)
    df = DataFrame(rows)
    sort!(df, [:pid])
    CSV.write(path, df)
    return path
end

function select_best_seed(summary_path::AbstractString)
    df = DataFrame(CSV.File(summary_path))
    isempty(df) && error("Cannot select a best seed from an empty PID summary.")
    sort!(df, [:median_identity, :mean_identity, :pid]; rev = [true, true, false])
    row = first(eachrow(df))
    return SeedSelection(;
        pid = Float64(row.pid),
        median_identity = Float64(row.median_identity),
        mean_identity = Float64(row.mean_identity),
        stockholm_path = String(row.stockholm_path),
        fasta_path = row.fasta_path === missing ? nothing : String(row.fasta_path),
        summary_path
    )
end

function build_thoraxe_msa(target::ResolvedTarget, workdir::AbstractString;
        pid_thresholds::AbstractVector{<:Real} = DEFAULT_PID_THRESHOLDS,
        specieslist::Union{Nothing, AbstractString} = nothing,
        cached_thoraxe_input_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        orthology::AbstractString = "1:1",
        specieslist_filter::Bool = true,
        biomart_datasets_filter::Bool = true,
        transcript_query_timeout_seconds::Union{Nothing, Real} = 180,
        transcript_query_timeout_max_seconds::Union{Nothing, Real} = 240,
        transcript_query_retries::Integer = 2,
        allow_specieslist_timeout_fallback::Bool = true,
        thoraxe_timeout_seconds::Union{Nothing, Real} = nothing)
    _orthology_relationships(orthology)
    # Species filters run only when transcript_query will create new input.
    if cached_thoraxe_input_dir === nothing && specieslist_filter
        species_filter = _resolve_effective_specieslist(target, specieslist, orthology)
    else
        species_filter = (specieslist = _normalized_specieslist(specieslist),
            warnings = String[])
    end
    if cached_thoraxe_input_dir === nothing && biomart_datasets_filter
        biomart_filter = _resolve_biomart_datasets_specieslist(
            target, species_filter.specieslist)
    else
        biomart_filter = (specieslist = species_filter.specieslist, warnings = String[])
    end
    effective_specieslist = biomart_filter.specieslist
    input_dir = _ensure_transcript_query(target, workdir;
        specieslist = effective_specieslist, overwrite,
        cached_input_dir = cached_thoraxe_input_dir,
        timeout_seconds = transcript_query_timeout_seconds,
        timeout_max_seconds = transcript_query_timeout_max_seconds,
        max_retries = transcript_query_retries,
        orthology,
        allow_specieslist_timeout_fallback)
    thoraxe_dir = _ensure_baseline_thoraxe(
        target, input_dir, workdir; specieslist = effective_specieslist, overwrite,
        timeout_seconds = thoraxe_timeout_seconds)
    msa,
    species = assemble_transcript_msa(thoraxe_dir, target.ensembl_gene_id, target.transcript_id)
    warnings = vcat(species_filter.warnings,
        biomart_filter.warnings,
        _biomart_transcript_query_warnings(input_dir, _thoraxe_logs_dir(workdir)),
        _validate_transcript_translation(target, msa))
    baseline_fasta, baseline_sto,
    sequence_fasta, species_file = _save_baseline_msa(workdir, msa, species; overwrite)

    summary_path = joinpath(_thoraxe_msa_dir(workdir), "best_seed.csv")
    score_rows = NamedTuple[]
    for pid in Float64.(pid_thresholds)
        fasta_path,
        sto_path = _generate_pid_seed(
            target, input_dir, workdir, pid, species_file;
            overwrite, timeout_seconds = thoraxe_timeout_seconds)
        identity = compute_identity_against_reference(baseline_fasta, fasta_path;
            logs_dir = joinpath(_thoraxe_logs_dir(workdir), "hhalign"),
            label = "pid$(format_pid(pid))")
        push!(score_rows,
            (;
                gene_id = target.ensembl_gene_id,
                transcript_id = target.transcript_id,
                pid,
                mean_identity = identity,
                median_identity = identity,
                n_samples = 1,
                n_sequences_msa0 = nsequences(msa),
                fasta_path,
                stockholm_path = sto_path
            ))
    end
    _summarize_pid_scores(score_rows, summary_path)
    best = select_best_seed(summary_path)

    return ThorAxeMSAResult(;
        input_dir,
        thoraxe_dir,
        msa_dir = _thoraxe_msa_dir(workdir),
        baseline_fasta,
        baseline_stockholm = baseline_sto,
        sequence_fasta,
        species_file,
        pid_summary = summary_path,
        best_seed = best,
        logs_dir = _thoraxe_logs_dir(workdir),
        warnings,
        status = :ok
    )
end

end
