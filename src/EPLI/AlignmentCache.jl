const _ALIGNMENT_CACHE_SCHEMA_VERSION = 4

function _aligner_provenance(aligner_cache_id::Union{Nothing, AbstractString};
        implementation_fn::Function = _implementation_identity)
    if aligner_cache_id === nothing
        @info "EPLI alignment reuse disabled: no aligner_cache_id was supplied."
        return nothing
    end
    isempty(strip(aligner_cache_id)) &&
        throw(ArgumentError("aligner_cache_id must be a nonblank string or nothing."))
    implementation = implementation_fn()
    implementation === nothing && return nothing
    # Normalize once so cached JSON and the in-memory record have the same types.
    return (; aligner_cache_id = String(aligner_cache_id),
        implementation = JSON.parse(JSON.json(implementation)))
end

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

function _alignment_cache_identity(input_fasta::AbstractString,
        run_label::AbstractString,
        aligner_args::Cmd,
        provenance)
    provenance === nothing && return nothing
    return Dict(
        "schema_version" => _ALIGNMENT_CACHE_SCHEMA_VERSION,
        "input_sha256" => _file_sha256(input_fasta),
        "aligner_cache_id" => provenance.aligner_cache_id,
        "implementation" => provenance.implementation,
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

function _alignment_output_matches(path::AbstractString, digest)
    try
        return _file_sha256(path) == digest
    catch err
        err isa SystemError || err isa Base.IOError || rethrow()
        return false
    end
end

function _cached_alignment_path(output_fasta::AbstractString, identity::AbstractDict)
    state = _read_alignment_cache(output_fasta)
    state isa AbstractDict || return nothing
    for key in ("schema_version", "input_sha256", "aligner_cache_id", "implementation",
        "aligner_args", "run_label")
        value = get(state, key, nothing)
        value isa typeof(identity[key]) && value == identity[key] || return nothing
    end
    get(state, "output_fasta", nothing) == basename(output_fasta) || return nothing
    recorded_path = get(state, "alignment_path", nothing)
    recorded_path isa AbstractString || return nothing
    occursin('\0', recorded_path) && return nothing
    alignment_path = normpath(joinpath(dirname(output_fasta), recorded_path))
    isfile(alignment_path) || return nothing
    _alignment_output_matches(alignment_path, get(state, "output_sha256", nothing)) ||
        return nothing
    return alignment_path
end

function _write_alignment_cache(output_fasta::AbstractString,
        alignment_path::AbstractString,
        identity::AbstractDict)
    state = copy(identity)
    state["alignment_path"] = relpath(abspath(alignment_path), dirname(abspath(output_fasta)))
    state["output_fasta"] = basename(output_fasta)
    state["output_sha256"] = _file_sha256(alignment_path)
    cache_path = _alignment_cache_path(output_fasta)
    tmp_path = string(cache_path, ".tmp")
    write_json(tmp_path, state)
    mv(tmp_path, cache_path; force = true)
    return nothing
end

function _alignment_fasta_records(path::AbstractString)
    # Read original rows: an MSA reader may pad rows or normalize unknown residues.
    parts = split(read(path, String), '>')
    isempty(strip(first(parts))) || error("FASTA sequence precedes its header: $(path).")
    records = Dict{String, String}()
    for part in Iterators.drop(parts, 1)
        fields = split(part, '\n'; limit = 2)
        name = strip(first(fields))
        isempty(name) && error("FASTA record has an empty ID: $(path).")
        haskey(records, name) &&
            error("FASTA record ID $(repr(name)) is repeated: $(path).")
        length(fields) == 2 || error("FASTA record $(repr(name)) has no sequence: $(path).")
        sequence = uppercase(filter(!isspace, fields[2]))
        isempty(sequence) && error("FASTA record $(repr(name)) has no sequence: $(path).")
        records[name] = sequence
    end
    isempty(records) && error("FASTA file has no sequences: $(path).")
    return records
end

function _validate_alignment_output(input_records::AbstractDict, path::AbstractString)
    records = _alignment_fasta_records(path)
    all(==(length(first(values(records)))), length.(values(records))) ||
        error("Aligner output FASTA rows have different lengths: $(path).")
    keys(records) == keys(input_records) ||
        error("Aligner output FASTA IDs do not match the input: $(path).")
    for (name, sequence) in records
        replace(sequence, '-' => "", '.' => "") ==
        replace(input_records[name], '-' => "", '.' => "") ||
            error("Aligner changed residues for FASTA record $(repr(name)): $(path).")
    end
    return nothing
end

function _run_aligner(aligner_fn::Function,
        input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString},
        run_label::AbstractString,
        overwrite::Bool,
        aligner_args::Cmd,
        aligner_provenance = nothing)
    identity = _alignment_cache_identity(input_fasta, run_label, aligner_args, aligner_provenance)
    if !overwrite && identity !== nothing
        cached_path = _cached_alignment_path(output_fasta, identity)
        cached_path === nothing || return cached_path
    end
    mkpath(dirname(output_fasta))
    # A failed writer must not leave its partial output certified by an older success.
    rm(_alignment_cache_path(output_fasta); force = true)
    input_records = _alignment_fasta_records(input_fasta)
    result = aligner_fn(input_fasta, output_fasta; logs_dir, run_label, aligner_args)
    path = _alignment_result_path(result, output_fasta)
    isfile(path) || error("Aligner did not write the expected MSA FASTA at $(path).")
    _validate_alignment_output(input_records, path)
    identity === nothing || _write_alignment_cache(output_fasta, path, identity)
    return path
end
