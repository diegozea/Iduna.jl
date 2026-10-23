# Stage Paths
# -----------

"""
    _pipeline_stage_dir(workdir, stage_key) -> String

Return the package-owned directory used to store state for one pipeline stage.

# Arguments

- `workdir::AbstractString`: Iduna work directory.
- `stage_key::AbstractString`: stable key identifying one stage instance.
"""
function _pipeline_stage_dir(workdir::AbstractString, stage_key::AbstractString)
    safe_key = replace(String(stage_key), ':' => "__", '/' => "__", '\\' => "__")
    return joinpath(workdir, ".iduna", "stages", safe_key)
end

"""
    _stage_state_path(stage_dir) -> String

Return the path to the JSON state file stored inside a stage directory.

# Arguments

- `stage_dir::AbstractString`: stage state directory.
"""
_stage_state_path(stage_dir::AbstractString) = joinpath(stage_dir, _STAGE_STATE_FILE)

_stage_output_exists(::Nothing) = true
_stage_output_exists(path::AbstractString) = ispath(path)
_stage_output_exists(paths::AbstractVector) = all(_stage_output_exists, paths)
_stage_output_exists(outputs::NamedTuple) = all(_stage_output_exists, values(outputs))

function _stage_output_exists(outputs::AbstractDict)
    return all(_stage_output_exists, values(outputs))
end

# Output Serialization
# --------------------

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

# Required files are separate from descriptive outputs: derived tables may be
# regenerated, and directories must never be recursively hashed as artifacts.
_stage_output_files(::Nothing) = String[]
_stage_output_files(path::AbstractString) = [String(path)]
function _stage_output_files(outputs::Union{NamedTuple, AbstractDict, AbstractVector})
    return reduce(vcat, (_stage_output_files(value) for value in values(outputs));
        init = String[])
end

_output_mapping_matches(_expected, _stored) = false
_output_mapping_matches(::Nothing, stored) = stored === nothing
_output_mapping_matches(expected::AbstractString, stored) = expected == stored

function _output_mapping_matches(expected::AbstractVector, stored)
    return stored isa AbstractVector && length(expected) == length(stored) &&
           all(_output_mapping_matches(a, b) for (a, b) in zip(expected, stored))
end

function _output_mapping_matches(expected::Union{NamedTuple, AbstractDict}, stored)
    return stored isa AbstractDict &&
           all(haskey(stored, String(key)) &&
               _output_mapping_matches(value, stored[String(key)])
    for (key, value) in pairs(expected))
end

function _stage_output_digests(outputs, workdir::AbstractString)
    paths = unique(_stage_output_files(outputs))
    return Dict(path => _file_sha256(_resolve_artifact_path(path, workdir))
    for path in paths)
end

function _stage_output_digest_matches(path::AbstractString, digest, workdir::AbstractString)
    digest isa AbstractString && occursin(r"^[0-9a-f]{64}$", digest) || return false
    file = _resolve_artifact_path(path, workdir)
    return isfile(file) && _file_sha256(file) == digest
end

function _stage_output_integrity_matches(state, stage_dir::AbstractString, requested)
    root = get(state, "output_root", nothing)
    root isa AbstractString && !isabspath(root) || return false
    workdir = normpath(joinpath(stage_dir, root))
    required = get(state, "required_outputs", nothing)
    digests = get(state, "output_digests", nothing)
    digests isa AbstractDict || return false
    _output_mapping_matches(required, get(state, "outputs", nothing)) || return false
    _output_mapping_matches(_relative_output_value(requested, workdir), required) ||
        return false
    paths = unique(_stage_output_files(required))
    Set(keys(digests)) == Set(paths) || return false
    return all(_stage_output_digest_matches(path, digests[path], workdir) for path in paths)
end

# State Reading
# -------------

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
    elseif !(state isa AbstractDict)
        return "state must be a JSON object"
    end
    return nothing
end

function _stage_existing_started_at(
        stage_dir::AbstractString, identity_hash::AbstractString)
    state = _read_stage_state(stage_dir)
    state isa AbstractDict || return nothing
    get(state, "identity_hash", nothing) == identity_hash || return nothing
    started = get(state, "started_at", nothing)
    return started isa AbstractString ? String(started) : nothing
end

# State Writing and Reuse
# -----------------------

"""
    _write_stage_state(stage_dir; stage, stage_key, status, identity,
                       outputs=NamedTuple(), required_outputs=outputs,
                       warnings=String[], exception=nothing,
                       action=nothing, workdir=stage_dir,
                       preserve_started_at=false, extra=NamedTuple())

Write the JSON state file that lets Iduna decide whether a stage can be reused.

This is an internal contract shared by pipeline stages. Paths under `workdir`
are stored as relative paths so a result directory can be moved.

# Arguments

- `stage_dir::AbstractString`: directory where the stage state file is written.

# Keywords

- `stage::AbstractString`: human-readable stage name.
- `stage_key::AbstractString`: stable key identifying this stage instance.
- `status::Symbol`: stage status, such as `:running`, `:done`, or `:failed`.
- `identity`: JSON-serializable data that defines the run inputs.
- `outputs = NamedTuple()`: JSON-serializable output paths or path groups.
- `required_outputs = outputs`: files whose contents must remain unchanged for
  reuse. Exclude directories, logs and regenerable derived artifacts.
- `warnings::AbstractVector{<:AbstractString} = String[]`: non-fatal warnings to
  record.
- `exception = nothing`: exception summary to record for failed stages.
- `action = nothing`: stage action, such as `:run`, `:rebuild`, or `:reuse`.
- `workdir::AbstractString = stage_dir`: root used to store output paths as
  relative paths.
- `preserve_started_at::Bool = false`: keep the previous start timestamp when
  updating a running stage with the same identity.
- `extra = NamedTuple()`: extra JSON-serializable fields to merge into the state.
"""
function _write_stage_state(stage_dir::AbstractString;
        stage::AbstractString,
        stage_key::AbstractString,
        status::Symbol,
        identity,
        outputs = NamedTuple(),
        required_outputs = outputs,
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
    relative_outputs = _relative_output_value(outputs, workdir)
    required = _relative_output_value(required_outputs, workdir)
    _output_mapping_matches(required, relative_outputs) ||
        error("Required stage files must match their declared output paths.")
    digests = status === :done ? _stage_output_digests(required, workdir) : nothing
    state = merge(
        (;
            schema_version = _STAGE_STATE_SCHEMA_VERSION,
            stage = String(stage),
            stage_key = String(stage_key),
            status = String(status),
            action = action === nothing ? nothing : String(action),
            identity,
            identity_hash = hash,
            outputs = relative_outputs,
            output_root = relpath(abspath(workdir), abspath(stage_dir)),
            required_outputs = required,
            output_digests = digests,
            warnings = String.(warnings),
            exception,
            started_at,
            updated_at = timestamp,
            finished_at
        ),
        extra)
    tmp_path = string(state_path, ".tmp")
    write_json(tmp_path, state)
    # Replace the state file only after the new JSON has been written completely.
    mv(tmp_path, state_path; force = true)
    return state_path
end

function _stage_manifest_mismatch(state, expected_hash, stage_label)
    if get(state, "schema_version", nothing) !== _STAGE_STATE_SCHEMA_VERSION ||
       !(get(state, "status", nothing) isa AbstractString)
        return (; reusable = false, status = :stale,
            warning = "$(stage_label) state has an unsupported schema or status; rebuilding.")
    end
    status = Symbol(state["status"])
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
    if get(state, "identity_hash", nothing) != expected_hash
        return (;
            reusable = false,
            status = :stale,
            warning = "$(stage_label) inputs changed; rebuilding stale cached outputs.")
    end
    return nothing
end

"""
    _classify_stage_state(stage_dir, identity, required_outputs; stage_label="stage")

Compare a saved stage state with the requested inputs and outputs.

# Arguments

- `stage_dir::AbstractString`: directory that may contain a stage state file.
- `identity`: JSON-serializable data that defines the requested run inputs.
- `required_outputs`: output paths or path groups that must exist for reuse.

# Keywords

- `stage_label::AbstractString = "stage"`: label used in warnings.

# Returns

- A named tuple describing whether the cached stage is reusable, missing, stale,
  incomplete, or unreadable.
"""
function _classify_stage_state(stage_dir::AbstractString,
        identity,
        required_outputs;
        stage_label::AbstractString = "stage")
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
    mismatch = _stage_manifest_mismatch(state, _identity_hash(identity), stage_label)
    mismatch === nothing || return mismatch
    if !outputs_ready
        return (;
            reusable = false,
            status = :unfinished,
            warning = "$(stage_label) outputs are incomplete despite a done $(_STAGE_STATE_FILE); rebuilding.")
    end
    intact = try
        _stage_output_integrity_matches(state, stage_dir, required_outputs)
    catch err
        err isa InterruptException && rethrow()
        false
    end
    if !intact
        return (; reusable = false, status = :stale,
            warning = "$(stage_label) output integrity could not be verified; rebuilding.")
    end
    return (; reusable = true, status = :done, warning = nothing)
end
