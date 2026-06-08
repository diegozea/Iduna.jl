module Utils

using BioAlignments: AffineGapScoreModel, BLOSUM62, GlobalAlignment, alignment,
                     count_deletions, count_insertions, count_mismatches, pairalign
using BioSequences: LongAA
import CodecZlib
using Dates: UTC, now
import HTTP
import JSON
import SHA
using MIToS.MSA: AbstractMultipleSequenceAlignment, getannotcolumn, getannotfile,
                 ncolumns, sequencenames, setannotcolumn!, setannotfile!
using Printf: @sprintf

include("Types.jl")

export DEFAULT_PID_THRESHOLDS,
       ExpansionResult,
       IdunaResult,
       ResolvedTarget,
       SeedSelection,
       ThorAxeMSAResult,
       ValidationResult,
       decode_body,
       ensure_mmseqs_db,
       fasta_sequence,
       format_fasta,
       id_kind,
       is_ensembl_transcript_id,
       is_uniprot_id,
       prepare_output_dir,
       protein_alignment_stats,
       resolve_sequence_name,
       result_summary,
       run_logged,
       safe_rm,
       sequence_name_variants,
       strip_ensembl_version,
       write_fasta,
       write_json,
       write_text,
       format_pid,
       format_pid_dir,
       S_EXON_CODE_FEATURE,
       S_EXON_CODE_MAP_FEATURE,
       S_EXON_MISSING_CODE,
       s_exon_blocks_path,
       s_exon_code_map,
       s_exon_codes,
       has_s_exon_annotations,
       set_s_exon_annotations!,
       write_s_exon_blocks_tsv

const DEFAULT_PID_THRESHOLDS = Float64[10, 20, 30, 60, 80]
const S_EXON_CODE_FEATURE = "SExonCode"
const S_EXON_CODE_MAP_FEATURE = "SExonCodeMap"
const S_EXON_MISSING_CODE = '.'
const _STAGE_STATE_SCHEMA_VERSION = 1
const _STAGE_STATE_FILE = "stage_state.json"
const _PROTEIN_ALIGNMENT_SCORE_MODEL = AffineGapScoreModel(BLOSUM62, gap_open = -10, gap_extend = -1)
const _TRANSIENT_HTTP_STATUSES = Set([429, 500, 502, 503, 504])

struct _RetryableHTTPStatus <: Exception
    response::HTTP.Response
end

# Julia's retry helper retries exceptions, so transient statuses are wrapped.
_is_transient_http_status(status::Integer)::Bool = status in _TRANSIENT_HTTP_STATUSES
_is_retryable_http_exception(_state, err)::Bool = err isa _RetryableHTTPStatus

_http_get_request(url::AbstractString; kwargs...) = HTTP.request("GET", url; kwargs...)

function _http_get_with_retries(url::AbstractString,
        headers;
        retries::Integer = 4,
        sleep_seconds::Real = 1.5,
        max_delay::Real = 30.0,
        http_get::Function = _http_get_request)
    attempts = max(Int(retries), 1)
    try
        # Return the final HTTP response so callers decide what status means success.
        return retry(;
            delays = Base.ExponentialBackOff(;
                n = attempts - 1,
                first_delay = Float64(sleep_seconds),
                max_delay = Float64(max_delay),
                factor = 2.0,
                jitter = 0.0),
            check = _is_retryable_http_exception) do
            resp = http_get(url; headers = headers, retry = false, status_exception = false)
            _is_transient_http_status(resp.status) && throw(_RetryableHTTPStatus(resp))
            return resp
        end()
    catch err
        err isa _RetryableHTTPStatus && return err.response
        rethrow()
    end
end

strip_ensembl_version(id::AbstractString)::String = String(split(String(id), '.'; limit = 2)[1])

function _strip_numeric_prefix(id::AbstractString)
    parts = split(String(id), '|'; limit = 2)
    return length(parts) == 2 && all(isdigit, parts[1]) ? String(parts[2]) : String(id)
end

function sequence_name_variants(id::AbstractString)
    variants = String[]
    for base in (String(id), _strip_numeric_prefix(id))
        push!(variants, base)
        core = strip_ensembl_version(base)
        core == base || push!(variants, core)
    end
    return unique(variants)
end

function resolve_sequence_name(msa::AbstractMultipleSequenceAlignment,
        ids;
        fallback::Bool = false)
    names = String.(sequencenames(msa))
    lookup = Dict{String, String}()
    # Store both full and version-stripped names for Ensembl IDs.
    for name in names
        for variant in sequence_name_variants(name)
            haskey(lookup, variant) || (lookup[variant] = name)
        end
    end

    iterable_ids = ids isa AbstractString ? (ids,) : ids
    for id in iterable_ids
        for variant in sequence_name_variants(id)
            haskey(lookup, variant) && return lookup[variant]
        end
    end
    return fallback && !isempty(names) ? first(names) : nothing
end

function protein_alignment_stats(query_seq::AbstractString,
        reference_seq::AbstractString;
        include_alignment::Bool = false)
    result = pairalign(GlobalAlignment(), LongAA(uppercase(String(reference_seq))),
        LongAA(uppercase(String(query_seq))), _PROTEIN_ALIGNMENT_SCORE_MODEL)
    aln = alignment(result)
    mismatches = count_mismatches(aln)
    insertions = count_insertions(aln)
    deletions = count_deletions(aln)
    stats = (;
        identical = mismatches == 0 && insertions == 0 && deletions == 0,
        mismatches,
        insertions,
        deletions
    )
    return include_alignment ? merge((; aln), stats) : stats
end

is_ensembl_transcript_id(id::AbstractString)::Bool = occursin(r"^ENS[A-Z]*T[0-9]+(\.[0-9]+)?$", String(id))

# UniProt documents accessions as either [OPQ][0-9][A-Z0-9]{3}[0-9] or
# [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}. The optional suffix accepts
# UniProt isoform IDs, which are written as accession-number plus "-<number>".
is_uniprot_id(id::AbstractString)::Bool = occursin(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9])(-[0-9]+)?$",
    String(id))

function id_kind(id::AbstractString)::Symbol
    is_ensembl_transcript_id(id) && return :ensembl_transcript
    is_uniprot_id(id) && return :uniprot
    error("Input ID $(id) is not recognized as a UniProt accession or Ensembl transcript ID.")
end

format_pid(pid::Real) = @sprintf("%.1f", Float64(pid))
format_pid_dir(pid::Real) = "pid_$(@sprintf("%.2f", Float64(pid)))"

function _body_bytes(body)::Vector{UInt8}
    body isa AbstractVector{UInt8} && return Vector{UInt8}(body)
    return Vector{UInt8}(codeunits(String(body)))
end

function decode_body(resp::HTTP.Response)::String
    body = _body_bytes(resp.body)
    if length(body) >= 2 && body[1] == 0x1f && body[2] == 0x8b
        return String(transcode(CodecZlib.GzipDecompressor, body))
    end
    return String(body)
end

function prepare_output_dir(input_id::AbstractString;
        workdir::Union{Nothing, AbstractString} = nothing,
        output_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false)
    if workdir !== nothing && output_dir !== nothing &&
       abspath(workdir) != abspath(output_dir)
        error("Use either workdir or output_dir, or pass the same path for both.")
    end
    root = workdir === nothing ? output_dir : workdir
    root = root === nothing ? abspath(String(input_id)) : abspath(String(root))

    if isfile(root)
        error("Output path $(root) exists and is a file.")
    end
    # The workdir itself is never deleted. When overwrite=true, each pipeline
    # stage removes only the package-owned subdirectories it is about to rebuild.
    mkpath(root)
    return root
end

function safe_rm(path::AbstractString, root::AbstractString)
    abs_path = abspath(path)
    abs_root = abspath(root)
    rel = relpath(abs_path, abs_root)
    # Deletions are allowed only inside the active work directory.
    if rel == "." || startswith(rel, "..") || isabspath(rel)
        error("Refusing to remove $(abs_path) because it is outside workdir $(abs_root).")
    end
    rm(abs_path; recursive = true, force = true)
    return nothing
end

function run_logged(cmd::Cmd;
        stdout_path::AbstractString,
        stderr_path::AbstractString,
        workdir::Union{Nothing, AbstractString} = nothing)
    mkpath(dirname(stdout_path))
    mkpath(dirname(stderr_path))
    open(stdout_path, "w") do out_io
        open(stderr_path, "w") do err_io
            if workdir === nothing
                run(pipeline(cmd; stdout = out_io, stderr = err_io))
            else
                cd(workdir) do
                    run(pipeline(cmd; stdout = out_io, stderr = err_io))
                end
            end
        end
    end
    return nothing
end

function fasta_sequence(content::AbstractString)::Union{Nothing, String}
    seq = String[]
    for line in split(content, '\n')
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, '>') && continue
        push!(seq, stripped)
    end
    isempty(seq) && return nothing
    return uppercase(join(seq))
end

function _wrap_sequence(seq::AbstractString; width::Int = 60)
    io = IOBuffer()
    i = firstindex(seq)
    while i <= lastindex(seq)
        j = min(i + width - 1, lastindex(seq))
        println(io, seq[i:j])
        i = j + 1
    end
    return String(take!(io))
end

format_fasta(id::AbstractString,
    seq::AbstractString)::String = string(">", id, "\n", _wrap_sequence(uppercase(String(seq))))

function write_fasta(path::AbstractString, records)
    mkpath(dirname(path))
    open(path, "w") do io
        for (name, seq) in records
            println(io, '>', name)
            println(io, _wrap_sequence(uppercase(String(seq))))
        end
    end
    return path
end

function write_text(path::AbstractString, text::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        print(io, text)
    end
    return path
end

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.json(io, obj; pretty = true)
        println(io)
    end
    return path
end

function _file_sha256(path::AbstractString)
    return bytes2hex(open(SHA.sha256, path))
end

function _identity_hash(identity)
    return bytes2hex(SHA.sha256(JSON.json(identity)))
end

function _pipeline_stage_dir(workdir::AbstractString, stage_key::AbstractString)
    safe_key = replace(String(stage_key), ':' => "__", '/' => "__", '\\' => "__")
    return joinpath(workdir, ".iduna", "stages", safe_key)
end

_stage_state_path(stage_dir::AbstractString) = joinpath(stage_dir, _STAGE_STATE_FILE)

_stage_output_exists(::Nothing) = true
_stage_output_exists(path::AbstractString) = ispath(path)
_stage_output_exists(paths::AbstractVector) = all(_stage_output_exists, paths)
_stage_output_exists(outputs::NamedTuple) = all(_stage_output_exists, values(outputs))

function _stage_output_exists(outputs::AbstractDict)
    return all(_stage_output_exists, values(outputs))
end

function _relative_output_value(path::AbstractString, workdir::AbstractString)
    return _relative_artifact_path(path, workdir)
end

_relative_output_value(::Nothing, _workdir::AbstractString) = nothing

function _relative_output_value(paths::AbstractVector, workdir::AbstractString)
    return [_relative_output_value(path, workdir) for path in paths]
end

function _relative_output_value(outputs::NamedTuple, workdir::AbstractString)
    return Dict(String(name) => _relative_output_value(value, workdir)
    for (name, value) in pairs(outputs))
end

function _relative_output_value(outputs::AbstractDict, workdir::AbstractString)
    return Dict(String(name) => _relative_output_value(value, workdir)
    for (name, value) in pairs(outputs))
end

function _read_stage_state(stage_dir::AbstractString)
    state_path = _stage_state_path(stage_dir)
    isfile(state_path) || return nothing
    try
        return JSON.parse(read(state_path, String))
    catch err
        err isa InterruptException && rethrow()
        return (; unreadable = sprint(showerror, err))
    end
end

function _stage_state_unreadable_message(state)
    if state === nothing
        return "state file disappeared while reading"
    elseif state isa NamedTuple && haskey(state, :unreadable)
        return state.unreadable
    end
    return nothing
end

function _stage_existing_started_at(stage_dir::AbstractString, identity_hash::AbstractString)
    state = _read_stage_state(stage_dir)
    state isa AbstractDict || return nothing
    String(get(state, "identity_hash", "")) == identity_hash || return nothing
    started = get(state, "started_at", nothing)
    return started isa AbstractString ? String(started) : nothing
end

function _write_stage_state(stage_dir::AbstractString;
        stage::AbstractString,
        stage_key::AbstractString,
        status::Symbol,
        identity,
        outputs = NamedTuple(),
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing,
        action::Union{Nothing, Symbol, AbstractString} = nothing,
        workdir::AbstractString = stage_dir,
        preserve_started_at::Bool = false,
        extra = NamedTuple())
    mkpath(stage_dir)
    state_path = _stage_state_path(stage_dir)
    hash = _identity_hash(identity)
    timestamp = string(now(UTC))
    started_at = _stage_existing_started_at(stage_dir, hash)
    if started_at === nothing || (status === :running && !preserve_started_at)
        started_at = timestamp
    end
    finished_at = status in (:done, :failed, :skipped) ? timestamp : nothing
    state = merge(
        (;
            schema_version = _STAGE_STATE_SCHEMA_VERSION,
            stage = String(stage),
            stage_key = String(stage_key),
            status = String(status),
            action = action === nothing ? nothing : String(action),
            identity,
            identity_hash = hash,
            outputs = _relative_output_value(outputs, workdir),
            warnings = String.(warnings),
            exception,
            started_at,
            updated_at = timestamp,
            finished_at
        ),
        extra)
    tmp_path = string(state_path, ".tmp")
    write_json(tmp_path, state)
    mv(tmp_path, state_path; force = true)
    return state_path
end

function _classify_stage_state(stage_dir::AbstractString,
        identity,
        required_outputs;
        stage_label::AbstractString = "stage")
    expected_hash = _identity_hash(identity)
    outputs_ready = _stage_output_exists(required_outputs)
    state_path = _stage_state_path(stage_dir)
    if !isfile(state_path)
        if outputs_ready
            return (;
                reusable = false,
                status = :stale,
                warning = "Existing $(stage_label) outputs have no $(_STAGE_STATE_FILE); rebuilding to verify run identity.")
        elseif isdir(stage_dir)
            return (;
                reusable = false,
                status = :unfinished,
                warning = "Previous $(stage_label) state directory has no $(_STAGE_STATE_FILE) and incomplete outputs; rebuilding.")
        end
        return (; reusable = false, status = :missing, warning = nothing)
    end

    state = _read_stage_state(stage_dir)
    unreadable = _stage_state_unreadable_message(state)
    if unreadable !== nothing
        return (;
            reusable = false,
            status = :stale,
            warning = "Could not read $(stage_label) $(_STAGE_STATE_FILE): $(unreadable); rebuilding.")
    end

    status = Symbol(String(get(state, "status", "stale")))
    if status == :running
        return (;
            reusable = false,
            status = :unfinished,
            warning = "Previous $(stage_label) status was running; rebuilding unfinished outputs.")
    elseif status != :done
        return (;
            reusable = false,
            status,
            warning = "Previous $(stage_label) status was $(status); rebuilding.")
    end
    if String(get(state, "identity_hash", "")) != expected_hash
        return (;
            reusable = false,
            status = :stale,
            warning = "$(stage_label) inputs changed; rebuilding stale cached outputs.")
    end
    if !outputs_ready
        return (;
            reusable = false,
            status = :unfinished,
            warning = "$(stage_label) outputs are incomplete despite a done $(_STAGE_STATE_FILE); rebuilding.")
    end
    return (; reusable = true, status = :done, warning = nothing)
end

function s_exon_blocks_path(stockholm_path::AbstractString)
    root, _ = splitext(String(stockholm_path))
    return string(root, "_s_exon_blocks.tsv")
end

function s_exon_codes(msa::AbstractMultipleSequenceAlignment;
        feature::AbstractString = S_EXON_CODE_FEATURE)
    codes = getannotcolumn(msa, String(feature), "")
    if isempty(codes)
        return repeat(string(S_EXON_MISSING_CODE), ncolumns(msa))
    end
    length(codes) == ncolumns(msa) ||
        error("$(feature) has $(length(codes)) characters, but the MSA has $(ncolumns(msa)) columns.")
    return String(codes)
end

function has_s_exon_annotations(msa::AbstractMultipleSequenceAlignment;
        feature::AbstractString = S_EXON_CODE_FEATURE)
    codes = getannotcolumn(msa, String(feature), "")
    return !isempty(codes) && length(codes) == ncolumns(msa)
end

function s_exon_code_map(msa::AbstractMultipleSequenceAlignment;
        feature::AbstractString = S_EXON_CODE_MAP_FEATURE)
    raw = getannotfile(msa, String(feature), "")
    map = Dict{Char, String}()
    isempty(raw) && return map
    pos = firstindex(raw)
    for match in eachmatch(r"(\"(?:\\.|[^\"])*\")=>(\"(?:\\.|[^\"])*\")", raw)
        match_start = match.offset
        separator = pos == match_start ? "" : raw[pos:prevind(raw, match_start)]
        all(==(','), separator) ||
            error("Invalid $(feature) entry $(repr(raw)); expected quoted code=>s_exon_id pairs.")
        code = only(String(Meta.parse(match.captures[1])))
        code == S_EXON_MISSING_CODE &&
            error("$(S_EXON_MISSING_CODE) is reserved for columns without s-exon provenance.")
        map[code] = String(Meta.parse(match.captures[2]))
        pos = nextind(raw, match.offset + ncodeunits(match.match) - 1)
    end
    pos > lastindex(raw) ||
        error("Invalid $(feature) entry $(repr(raw)); expected quoted code=>s_exon_id pairs.")
    return map
end

function _format_s_exon_code_map(code_map)
    entries = ["$(repr(string(code)))=>$(repr(String(id)))" for (code, id) in code_map]
    return join(entries, ',')
end

function set_s_exon_annotations!(msa::AbstractMultipleSequenceAlignment,
        codes::AbstractString,
        code_map;
        code_feature::AbstractString = S_EXON_CODE_FEATURE,
        map_feature::AbstractString = S_EXON_CODE_MAP_FEATURE)
    length(codes) == ncolumns(msa) ||
        error("$(code_feature) has $(length(codes)) characters, but the MSA has $(ncolumns(msa)) columns.")
    setannotcolumn!(msa, String(code_feature), String(codes))
    setannotfile!(msa, String(map_feature), _format_s_exon_code_map(code_map))
    return msa
end

function _write_s_exon_block_line(io,
        alignment::AbstractString,
        pid_value::AbstractString,
        code::Char,
        s_exon_id::AbstractString,
        start_col::Integer,
        end_col::Integer)
    println(io,
        alignment, '\t',
        pid_value, '\t',
        code, '\t',
        s_exon_id, '\t',
        start_col, '\t',
        end_col, '\t',
        end_col - start_col + 1)
    return nothing
end

function _write_s_exon_block_rows(io,
        codes::AbstractString,
        code_map,
        alignment::AbstractString,
        pid_value::AbstractString)
    start_col = 0
    current = S_EXON_MISSING_CODE
    for (idx, code) in enumerate(codes)
        if code != current
            if current != S_EXON_MISSING_CODE
                _write_s_exon_block_line(io, alignment, pid_value, current,
                    get(code_map, current, ""), start_col, idx - 1)
            end
            current = code
            start_col = idx
        end
    end
    current == S_EXON_MISSING_CODE && return nothing
    _write_s_exon_block_line(io, alignment, pid_value, current,
        get(code_map, current, ""), start_col, length(codes))
    return nothing
end

function write_s_exon_blocks_tsv(path::AbstractString,
        msa::AbstractMultipleSequenceAlignment;
        alignment::AbstractString,
        pid = nothing,
        append::Bool = false,
        code_feature::AbstractString = S_EXON_CODE_FEATURE,
        map_feature::AbstractString = S_EXON_CODE_MAP_FEATURE)
    codes = s_exon_codes(msa; feature = code_feature)
    code_map = s_exon_code_map(msa; feature = map_feature)
    mkpath(dirname(path))
    pid_value = pid === nothing ? "" : string(pid)
    open(path, append ? "a" : "w") do io
        append ||
            println(io, "alignment\tpid\tcode\ts_exon_id\tstart_col\tend_col\tn_columns")
        _write_s_exon_block_rows(io, codes, code_map, alignment, pid_value)
    end
    return path
end

function ensure_mmseqs_db(db::AbstractString)
    for db_path in (String(db), string(db, "_aln"), string(db, "_seq"))
        dbtype = string(db_path, ".dbtype")
        isfile(dbtype) ||
            error("MMseqs2 database at $(db_path) looks incomplete; missing $(dbtype).")
    end
    return String(db)
end

function _is_path_inside_workdir(rel::AbstractString)
    rel == "." && return true
    isabspath(rel) && return false
    parts = splitpath(rel)
    return !isempty(parts) && first(parts) != ".."
end

_relative_artifact_path(path::Nothing, _workdir::AbstractString) = nothing

function _relative_artifact_path(path::AbstractString, workdir::AbstractString)
    str = String(path)
    isabspath(str) || return str
    rel = relpath(abspath(str), abspath(workdir))
    return _is_path_inside_workdir(rel) ? rel : str
end

_resolve_artifact_path(path::Nothing, _workdir::AbstractString) = nothing

function _resolve_artifact_path(path::AbstractString, workdir::AbstractString)
    str = String(path)
    return isabspath(str) ? str : joinpath(workdir, str)
end

function _resolve_artifact_path(result::IdunaResult, path)
    _resolve_artifact_path(path, result.workdir)
end

function _relative_seed_paths(seed::SeedSelection, workdir::AbstractString)
    return SeedSelection(;
        pid = seed.pid,
        median_identity = seed.median_identity,
        mean_identity = seed.mean_identity,
        stockholm_path = _relative_artifact_path(seed.stockholm_path, workdir),
        fasta_path = _relative_artifact_path(seed.fasta_path, workdir),
        s_exon_blocks_tsv = _relative_artifact_path(seed.s_exon_blocks_tsv, workdir),
        summary_path = _relative_artifact_path(seed.summary_path, workdir),
        used_fallback_dir = seed.used_fallback_dir,
        workdir
    )
end

function _relative_target_paths(target::ResolvedTarget, workdir::AbstractString)
    return ResolvedTarget(;
        input_id = target.input_id,
        input_kind = target.input_kind,
        uniprot_id = target.uniprot_id,
        ensembl_gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        ensembl_protein_id = target.ensembl_protein_id,
        species = target.species,
        uniprot_sequence_path = _relative_artifact_path(target.uniprot_sequence_path, workdir),
        ensembl_protein_sequence_path = _relative_artifact_path(
            target.ensembl_protein_sequence_path, workdir),
        sequence_validated = target.sequence_validated,
        mapping_confirmed = target.mapping_confirmed,
        workdir,
        warnings = target.warnings
    )
end

function _relative_thoraxe_msa_paths(thoraxe::ThorAxeMSAResult, workdir::AbstractString)
    return ThorAxeMSAResult(;
        input_dir = _relative_artifact_path(thoraxe.input_dir, workdir),
        thoraxe_dirs = [_relative_artifact_path(path, workdir)
                        for path in thoraxe.thoraxe_dirs],
        msa_dir = _relative_artifact_path(thoraxe.msa_dir, workdir),
        baseline_fastas = [_relative_artifact_path(path, workdir)
                           for path in thoraxe.baseline_fastas],
        baseline_stockholms = [_relative_artifact_path(path, workdir)
                               for path in thoraxe.baseline_stockholms],
        sequence_fastas = [_relative_artifact_path(path, workdir)
                           for path in thoraxe.sequence_fastas],
        species_files = [_relative_artifact_path(path, workdir)
                         for path in thoraxe.species_files],
        pid_summary = _relative_artifact_path(thoraxe.pid_summary, workdir),
        seeds = [_relative_seed_paths(seed, workdir) for seed in thoraxe.seeds],
        logs_dir = _relative_artifact_path(thoraxe.logs_dir, workdir),
        pid_sample_count = thoraxe.pid_sample_count,
        pid_sample_fraction = thoraxe.pid_sample_fraction,
        pid_sample_seed = thoraxe.pid_sample_seed,
        sampling_strategy = thoraxe.sampling_strategy,
        warnings = thoraxe.warnings,
        status = thoraxe.status
    )
end

function _relative_expansion_paths(expansion::ExpansionResult, workdir::AbstractString)
    return ExpansionResult(;
        run_dir = _relative_artifact_path(expansion.run_dir, workdir),
        seed_stockholm = _relative_artifact_path(expansion.seed_stockholm, workdir),
        seed_fasta = _relative_artifact_path(expansion.seed_fasta, workdir),
        hits_fasta = _relative_artifact_path(expansion.hits_fasta, workdir),
        full_stockholm = _relative_artifact_path(expansion.full_stockholm, workdir),
        match_stockholm = _relative_artifact_path(expansion.match_stockholm, workdir),
        a3m_path = _relative_artifact_path(expansion.a3m_path, workdir),
        s_exon_blocks_tsv = _relative_artifact_path(expansion.s_exon_blocks_tsv, workdir),
        db_dir = _relative_artifact_path(expansion.db_dir, workdir),
        hmm_dir = _relative_artifact_path(expansion.hmm_dir, workdir),
        logs_dir = _relative_artifact_path(expansion.logs_dir, workdir),
        n_hits = expansion.n_hits,
        n_new_hits = expansion.n_new_hits,
        status = expansion.status,
        workdir
    )
end

_relative_expansion_paths(::Missing, _workdir::AbstractString) = missing

function _relative_validation_paths(validation::ValidationResult, workdir::AbstractString)
    return ValidationResult(;
        stats_path = _relative_artifact_path(validation.stats_path, workdir),
        query_name = validation.query_name,
        query_vs_uniprot_path = _relative_artifact_path(
            validation.query_vs_uniprot_path, workdir),
        seed_nseq = validation.seed_nseq,
        seed_ncol = validation.seed_ncol,
        seed_clusters62 = validation.seed_clusters62,
        seed_neff80 = validation.seed_neff80,
        expanded_nseq = validation.expanded_nseq,
        expanded_ncol = validation.expanded_ncol,
        expanded_clusters62 = validation.expanded_clusters62,
        expanded_neff80 = validation.expanded_neff80,
        aln_identical = validation.aln_identical,
        aln_mismatches = validation.aln_mismatches,
        aln_insertions = validation.aln_insertions,
        aln_deletions = validation.aln_deletions,
        warnings = validation.warnings,
        status = validation.status
    )
end

function _relative_result_paths(result::IdunaResult)
    return IdunaResult(;
        input_id = result.input_id,
        workdir = result.workdir,
        target = _relative_target_paths(result.target, result.workdir),
        thoraxe_msa = _relative_thoraxe_msa_paths(result.thoraxe_msa, result.workdir),
        expansions = [_relative_expansion_paths(expansion, result.workdir)
                      for expansion in result.expansions],
        validations = [_relative_validation_paths(validation, result.workdir)
                       for validation in result.validations],
        stages = result.stages,
        warnings = result.warnings,
        status = result.status
    )
end

function _stage_summary_from_state_path(state_path::AbstractString, workdir::AbstractString)
    state = try
        JSON.parse(read(state_path, String))
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    state isa AbstractDict || return nothing
    return (;
        stage = get(state, "stage", nothing),
        stage_key = get(state, "stage_key", nothing),
        status = get(state, "status", nothing),
        action = get(state, "action", nothing),
        identity_hash = get(state, "identity_hash", nothing),
        state_path = _relative_artifact_path(state_path, workdir),
        outputs = get(state, "outputs", Dict{String, Any}()),
        warnings = get(state, "warnings", String[]),
        exception = get(state, "exception", nothing)
    )
end

function _push_stage_state_paths!(paths::Vector{String}, root::AbstractString)
    isdir(root) || return paths
    for (dir, _, files) in walkdir(root)
        _STAGE_STATE_FILE in files || continue
        push!(paths, joinpath(dir, _STAGE_STATE_FILE))
    end
    return paths
end

function collect_stage_summaries(workdir::AbstractString; stage_keys = nothing)
    paths = String[]
    _push_stage_state_paths!(paths, joinpath(workdir, ".iduna", "stages"))
    _push_stage_state_paths!(paths, joinpath(workdir, "expansion"))
    _push_stage_state_paths!(paths, joinpath(workdir, "validation"))
    stage_key_set = stage_keys === nothing ? nothing : Set(String.(stage_keys))
    summaries = Any[]
    for path in sort!(unique(paths))
        summary = _stage_summary_from_state_path(path, workdir)
        summary === nothing && continue
        if stage_key_set !== nothing
            summary.stage_key === nothing && continue
            String(summary.stage_key) in stage_key_set || continue
        end
        push!(summaries, summary)
    end
    return summaries
end

function _seed_summary(seed::SeedSelection)
    return (;
        pid = seed.pid,
        median_identity = seed.median_identity === missing ? nothing :
                          seed.median_identity,
        mean_identity = seed.mean_identity === missing ? nothing : seed.mean_identity,
        stockholm_path = seed.stockholm_path,
        fasta_path = seed.fasta_path,
        s_exon_blocks_tsv = seed.s_exon_blocks_tsv,
        summary_path = seed.summary_path,
        used_fallback_dir = seed.used_fallback_dir
    )
end

function _expansion_summary(expansion::ExpansionResult)
    return (;
        run_dir = expansion.run_dir,
        seed_stockholm = expansion.seed_stockholm,
        seed_fasta = expansion.seed_fasta,
        match_stockholm = expansion.match_stockholm,
        full_stockholm = expansion.full_stockholm,
        a3m_path = expansion.a3m_path,
        s_exon_blocks_tsv = expansion.s_exon_blocks_tsv,
        hits_fasta = expansion.hits_fasta,
        db_dir = expansion.db_dir,
        hmm_dir = expansion.hmm_dir,
        logs_dir = expansion.logs_dir,
        n_hits = expansion.n_hits,
        n_new_hits = expansion.n_new_hits,
        status = String(expansion.status)
    )
end

_expansion_summary(::Missing) = nothing

function _validation_summary(validation::ValidationResult)
    return (;
        stats_path = validation.stats_path,
        query_name = validation.query_name,
        query_vs_uniprot_path = validation.query_vs_uniprot_path,
        seed_nseq = validation.seed_nseq,
        seed_ncol = validation.seed_ncol,
        seed_clusters62 = validation.seed_clusters62,
        seed_neff80 = validation.seed_neff80,
        expanded_nseq = validation.expanded_nseq,
        expanded_ncol = validation.expanded_ncol,
        expanded_clusters62 = validation.expanded_clusters62,
        expanded_neff80 = validation.expanded_neff80,
        aln_identical = validation.aln_identical,
        aln_mismatches = validation.aln_mismatches,
        aln_insertions = validation.aln_insertions,
        aln_deletions = validation.aln_deletions,
        warnings = validation.warnings,
        status = String(validation.status)
    )
end

function result_summary(result::IdunaResult)
    result = _relative_result_paths(result)
    return (;
        input_id = result.input_id,
        status = String(result.status),
        warnings = result.warnings,
        stages = result.stages,
        target = (;
            input_kind = String(result.target.input_kind),
            uniprot_id = result.target.uniprot_id,
            ensembl_gene_id = result.target.ensembl_gene_id,
            transcript_id = result.target.transcript_id,
            ensembl_protein_id = result.target.ensembl_protein_id,
            species = result.target.species,
            uniprot_sequence_path = result.target.uniprot_sequence_path,
            ensembl_protein_sequence_path = result.target.ensembl_protein_sequence_path,
            sequence_validated = result.target.sequence_validated,
            mapping_confirmed = result.target.mapping_confirmed
        ),
        thoraxe_msa = (;
            baseline_fastas = result.thoraxe_msa.baseline_fastas,
            baseline_stockholms = result.thoraxe_msa.baseline_stockholms,
            pid_summary = result.thoraxe_msa.pid_summary,
            seeds = _seed_summary.(result.thoraxe_msa.seeds),
            selected_pids = [seed.pid for seed in result.thoraxe_msa.seeds],
            seed_stockholms = [seed.stockholm_path for seed in result.thoraxe_msa.seeds],
            s_exon_blocks_tsvs = [seed.s_exon_blocks_tsv
                                  for seed in result.thoraxe_msa.seeds],
            pid_sample_count = result.thoraxe_msa.pid_sample_count,
            pid_sample_fraction = result.thoraxe_msa.pid_sample_fraction,
            pid_sample_seed = result.thoraxe_msa.pid_sample_seed,
            sampling_strategy = String(result.thoraxe_msa.sampling_strategy),
            warnings = result.thoraxe_msa.warnings,
            status = String(result.thoraxe_msa.status)
        ),
        expansions = _expansion_summary.(result.expansions),
        validations = _validation_summary.(result.validations)
    )
end

end
