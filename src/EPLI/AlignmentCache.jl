const _ALIGNMENT_CACHE_SCHEMA_VERSION = 2

_aligner_cache_label(fn) = string(typeof(fn), "|", hash(fn))

function _alignment_result_path(result, fallback::AbstractString)
    result === nothing && return fallback
    result isa AbstractString && return String(result)
    if result isa NamedTuple && :fasta_path in keys(result)
        return String(result.fasta_path)
    end
    error("Aligner functions must return nothing, an output path, or a named tuple with `fasta_path`.")
end

function _alignment_cache_path(output_fasta::AbstractString)
    string(output_fasta, ".iduna_cache.json")
end

function _alignment_cache_identity(aligner_fn::Function,
        input_fasta::AbstractString,
        run_label::AbstractString,
        aligner_args::Cmd)
    return Dict(
        "schema_version" => _ALIGNMENT_CACHE_SCHEMA_VERSION,
        "input_sha256" => _file_sha256(input_fasta),
        "aligner" => _aligner_cache_label(aligner_fn),
        "aligner_args" => string(aligner_args),
        "run_label" => String(run_label))
end

function _read_alignment_cache(output_fasta::AbstractString)
    cache_path = _alignment_cache_path(output_fasta)
    isfile(cache_path) || return nothing
    try
        return JSON.parse(read(cache_path, String))
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _cached_alignment_path(output_fasta::AbstractString, identity::AbstractDict)
    state = _read_alignment_cache(output_fasta)
    state isa AbstractDict || return nothing
    for key in ("schema_version", "input_sha256", "aligner", "aligner_args", "run_label")
        get(state, key, nothing) == identity[key] || return nothing
    end
    alignment_path = get(state, "alignment_path", output_fasta)
    alignment_path isa AbstractString && isfile(alignment_path) || return nothing
    return String(alignment_path)
end

function _write_alignment_cache(output_fasta::AbstractString,
        alignment_path::AbstractString,
        identity::AbstractDict)
    state = copy(identity)
    state["alignment_path"] = String(alignment_path)
    state["output_fasta"] = String(output_fasta)
    write_json(_alignment_cache_path(output_fasta), state)
    return nothing
end

function _run_aligner(aligner_fn::Function,
        input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString},
        run_label::AbstractString,
        overwrite::Bool,
        aligner_args::Cmd)
    identity = _alignment_cache_identity(aligner_fn, input_fasta, run_label, aligner_args)
    if !overwrite
        cached_path = _cached_alignment_path(output_fasta, identity)
        cached_path === nothing || return cached_path
    end
    mkpath(dirname(output_fasta))
    result = aligner_fn(input_fasta, output_fasta; logs_dir, run_label, aligner_args)
    path = _alignment_result_path(result, output_fasta)
    isfile(path) || error("Aligner did not write the expected MSA FASTA at $(path).")
    _write_alignment_cache(output_fasta, path, identity)
    return path
end
