# Input Stage State
# -----------------

function _thoraxe_input_outputs(input_dir::AbstractString)
    ensembl_dir = joinpath(input_dir, "Ensembl")
    outputs = Dict{String, Any}()
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
    if status === :done
        _has_valid_ensembl_bundle(input_dir) ||
            error("Cannot complete an invalid Ensembl bundle.")
        metadata = JSON.parse(read(_transcript_query_metadata_path(input_dir), String))
        metadata isa AbstractDict || error("Transcript-query metadata must be an object.")
    end
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
    known = get(identity, :implementation_hash, nothing) !== nothing &&
            get(identity, :specieslist_fingerprint, nothing) !== nothing &&
            (get(identity, :source_kind, nothing) == "cached_input" ||
             get(identity, :tool_hash, nothing) !== nothing)
    known || return (; reusable = false, status = :stale,
        warning = "ThorAxe input producer identity is unknown; cache reuse is disabled.")
    return _classify_stage_state(_thoraxe_input_stage_dir(workdir), identity,
        _thoraxe_input_outputs(input_dir); stage_label = "ThorAxe transcript_query input")
end

function _maybe_reuse_thoraxe_input(workdir::AbstractString,
        input_dir::AbstractString,
        identity,
        cache;
        overwrite::Bool,
        manifest_message::AbstractString)
    overwrite && return nothing
    if cache.reusable
        @info manifest_message input_dir
        _write_thoraxe_input_state(workdir, input_dir, :done, identity; action = :reuse)
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
        source_fingerprint::Union{Nothing, AbstractString} = nothing,
        implementation = _thoraxe_implementation_identity(),
        tool_identity = nothing)
    return (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        species = _normalize_species_name(target.species),
        specieslist = _normalized_specieslist(specieslist),
        specieslist_fingerprint = _specieslist_identity_hash(specieslist),
        orthology = String(orthology),
        source_kind = String(source_kind),
        source_path = source_path === nothing ? nothing : abspath(String(source_path)),
        source_fingerprint,
        implementation_hash = implementation === nothing ? nothing :
                              _identity_hash(implementation),
        tool_hash = tool_identity === nothing ? nothing : _identity_hash(tool_identity)
    )
end

function _specieslist_identity_hash(specieslist::Union{Nothing, AbstractString})
    value = _normalized_specieslist(specieslist)
    return _inspect_thoraxe_identity() do
        content = value !== nothing && isfile(value) ? _file_sha256(value) : nothing
        _identity_hash((; value, content))
    end
end

function _thoraxe_reference_identity(target::ResolvedTarget, workdir::AbstractString)
    id = target.uniprot_id
    path = target.uniprot_sequence_path
    path === nothing && return (; id, status = :unrequested, sha256 = nothing)
    artifact_workdir = target.workdir === nothing ? workdir : target.workdir
    resolved = _resolve_artifact_path(path, artifact_workdir)
    isfile(resolved) || return (; id, status = :missing, sha256 = nothing)
    return _inspect_thoraxe_identity() do
        (; id, status = :present, sha256 = _file_sha256(resolved))
    end
end

function _metadata_value_matches(stored, expected)
    if expected === nothing
        return stored === nothing || stored === missing
    end
    stored === nothing && return false
    stored === missing && return false
    return String(stored) == String(expected)
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

function _missing_transcript_query_outputs(tmp_gene_dir::AbstractString)
    ensembl_dir = joinpath(tmp_gene_dir, "Ensembl")
    missing = String[]
    if isdir(ensembl_dir)
        for file in _REQUIRED_ENSEMBL_FILES
            path = joinpath(ensembl_dir, file)
            if !(isfile(path) && filesize(path) > 0)
                push!(missing, file)
            end
        end
    else
        append!(missing, _REQUIRED_ENSEMBL_FILES)
    end
    return missing
end
