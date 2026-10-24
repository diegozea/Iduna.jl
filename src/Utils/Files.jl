# Formatting and Response Decoding
# --------------------------------

"""
    format_pid(pid) -> String

Format a percent identity value for user-facing tables and messages.

# Arguments

- `pid::Real`: percent identity value.
"""
format_pid(pid::Real) = @sprintf("%.1f", Float64(pid))

"""
    format_pid_dir(pid) -> String

Format a percent identity threshold as a distinct directory name for its effective
`Float64` value. Keep two decimal places when they represent that value exactly;
otherwise use a spelling that preserves its full precision. Signed zero uses `pid_0.00`.

# Arguments

- `pid::Real`: finite percent identity threshold in `[0, 100]` after `Float64` conversion.

# Throws

- `ArgumentError`: if the effective threshold is non-finite or outside `[0, 100]`.
"""
function format_pid_dir(pid::Real)
    value = _effective_pid(pid)
    short = @sprintf("%.2f", value)
    token = parse(Float64, short) == value ? short : repr(value)
    return "pid_$(token)"
end

function _effective_pid(pid::Real)
    value = Float64(pid)
    isfinite(value) && 0.0 <= value <= 100.0 ||
        throw(ArgumentError("PID thresholds must be finite percentages in [0, 100]; got $(pid)."))
    return iszero(value) ? 0.0 : value
end

function _validate_pid_thresholds(pid_thresholds::AbstractVector{<:Real})
    isempty(pid_thresholds) && throw(ArgumentError("pid_thresholds cannot be empty."))
    values = _effective_pid.(pid_thresholds)
    allunique(values) ||
        throw(ArgumentError("pid_thresholds must be distinct after Float64 conversion."))
    return values
end

function _body_bytes(body)::Vector{UInt8}
    body isa AbstractVector{UInt8} && return Vector{UInt8}(body)
    return Vector{UInt8}(codeunits(String(body)))
end

"""
    decode_body(resp) -> String

Read an HTTP response body as text, including gzipped response bodies.

# Arguments

- `resp::HTTP.Response`: HTTP response whose body should be decoded.
"""
function decode_body(resp::HTTP.Response)::String
    body = _body_bytes(resp.body)
    if length(body) >= 2 && body[1] == 0x1f && body[2] == 0x8b
        return String(transcode(CodecZlib.GzipDecompressor, body))
    end
    return String(body)
end

# Workdir and Command Helpers
# ---------------------------

"""
    prepare_output_dir(input_id; workdir=nothing, output_dir=nothing, overwrite=false)

Choose and create the Iduna work directory.

# Arguments

- `input_id::AbstractString`: input ID used to derive the default output
  directory when no directory is passed.

# Keywords

- `workdir = nothing`: preferred name for the output directory. When `nothing`,
  `output_dir` is used, or the input ID is used if `output_dir` is also
  `nothing`.
- `output_dir = nothing`: older name for the same setting. When `nothing`,
  `workdir` controls the output directory.
- `overwrite::Bool = false`: accepted for caller consistency; stage cleanup is
  handled by each stage.

# Returns

- The absolute path to the work directory.
"""
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

function _is_within_directory(path::AbstractString, root::AbstractString)
    relative = relpath(path, root)
    return !isabspath(relative) && first(splitpath(relative)) != ".."
end

"""
    safe_rm(path, root)

Remove a path only when it and its resolved parent are inside the work directory.
The root itself is never removed. A symlink at `path` is removed without following
its target; parent symlinks must stay within the resolved root. Missing paths are
ignored only after checking their existing ancestors. The root may be a symlink.
This check assumes directories and symlinks are not changed concurrently.

# Arguments

- `path::AbstractString`: file or directory to remove.
- `root::AbstractString`: work directory that must contain `path`.

# Throws

- `ErrorException`: if `path` is the root, escapes it, or has an unresolved parent symlink.
"""
function safe_rm(path::AbstractString, root::AbstractString)
    # Strip trailing separators so rm unlinks a final directory symlink.
    abs_path = joinpath(splitpath(abspath(path)))
    abs_root = abspath(root)
    if relpath(abs_path, abs_root) == "." || !_is_within_directory(abs_path, abs_root)
        error("Refusing to remove $(abs_path): path must be strictly inside workdir $(abs_root).")
    end
    resolved_root = realpath(abs_root)
    # Resolve the parent, so removing a symlink itself does not follow its target.
    parent = dirname(abs_path)
    while !ispath(parent)
        islink(parent) &&
            error("Refusing to remove $(abs_path): ancestor $(parent) is an unresolved symlink.")
        parent = dirname(parent)
    end
    if !_is_within_directory(realpath(parent), resolved_root)
        error("Refusing to remove $(abs_path): parent resolves outside workdir $(abs_root).")
    end
    rm(abs_path; recursive = true, force = true)
    return nothing
end

"""
    run_logged(cmd; stdout_path, stderr_path, workdir=nothing)

Run an external command and write its standard output and standard error to
files.

# Arguments

- `cmd::Cmd`: command to run.

# Keywords

- `stdout_path::AbstractString`: file for standard output.
- `stderr_path::AbstractString`: file for standard error.
- `workdir = nothing`: directory where the command should run. When `nothing`,
  the command runs in the current directory.
"""
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

# FASTA and File Writers
# ----------------------

"""
    fasta_sequence(content) -> Union{Nothing, String}

Extract the sequence from FASTA text and return it in upper case.

# Arguments

- `content::AbstractString`: FASTA text to parse.
"""
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

"""
    format_fasta(id, seq) -> String

Format one protein sequence as FASTA text.

# Arguments

- `id::AbstractString`: FASTA record identifier.
- `seq::AbstractString`: protein sequence to write.
"""
format_fasta(id::AbstractString,
    seq::AbstractString)::String = string(">", id, "\n", _wrap_sequence(uppercase(String(seq))))

"""
    write_fasta(path, records) -> String

Write named sequences to a FASTA file.

# Arguments

- `path::AbstractString`: output file path.
- `records`: iterable of `(name, sequence)` pairs.
"""
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

"""
    write_text(path, text) -> String

Write text to a file and create parent directories when needed.

# Arguments

- `path::AbstractString`: output file path.
- `text::AbstractString`: text to write.
"""
function write_text(path::AbstractString, text::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        print(io, text)
    end
    return path
end

"""
    write_json(path, obj) -> String

Write an object as pretty JSON and create parent directories when needed.

# Arguments

- `path::AbstractString`: output file path.
- `obj`: JSON-serializable object to write.
"""
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

# Implementation Identity
# -----------------------

# Source and loaded dependency versions must stay stable during this Julia process.
# External executable bytes are checked separately for each public operation.
const _IMPLEMENTATION_IDENTITY_LOCK = ReentrantLock()
const _IMPLEMENTATION_SOURCE = Ref{Union{Nothing, String}}(nothing)

function _source_fingerprint(root::AbstractString)
    try
        files = String["Project.toml"]
        for directory in ("src", "ext")
            path = joinpath(root, directory)
            isdir(path) || continue
            for (parent, _, names) in walkdir(path), name in names

                # walkdir lists untraversed directory symlinks among file entries.
                isfile(joinpath(parent, name)) || return nothing
                occursin(r"\.jl(?:\.\d+)?\.cov$", name) && continue
                push!(files, relpath(joinpath(parent, name), root))
            end
        end
        records = [(path = file, sha256 = _file_sha256(joinpath(root, file)))
                   for file in sort!(files)]
        return _identity_hash(records)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _loaded_source_identity()
    return lock(_IMPLEMENTATION_IDENTITY_LOCK) do
        if _IMPLEMENTATION_SOURCE[] === nothing
            _IMPLEMENTATION_SOURCE[] = _source_fingerprint(dirname(dirname(@__DIR__)))
        end
        _IMPLEMENTATION_SOURCE[]
    end
end

# Upgraded or path-based stdlibs need their own package source identity.
function _is_bundled_stdlib(uuid, info)
    return Pkg.Types.is_stdlib(uuid) && !info.is_tracking_path &&
           info.tree_hash === nothing
end

function _dependency_record(uuid, info)
    stdlib = _is_bundled_stdlib(uuid, info)
    source_sha256 = nothing
    if !stdlib
        info.version === nothing && return nothing
        if info.is_tracking_path
            source_sha256 = _source_fingerprint(info.source)
            source_sha256 === nothing && return nothing
        else
            info.tree_hash === nothing && return nothing
            isdir(info.source) || return nothing
        end
    end
    return (;
        uuid = string(uuid),
        name = info.name,
        version = info.version === nothing ? nothing : string(info.version),
        tree_hash = info.tree_hash,
        source_sha256,
        julia = stdlib ? string(VERSION) : nothing
    )
end

function _dependency_module_matches(mod, uuid, info, record)
    Base.PkgId(mod) == Base.PkgId(uuid, info.name) || return false
    Base.pkgversion(mod) == info.version || return false
    source = pkgdir(mod)
    source === nothing && return record.julia !== nothing
    if isdir(source) && isdir(info.source)
        return realpath(source) == realpath(info.source)
    end
    return record.julia !== nothing && normpath(source) == normpath(info.source)
end

function _dependency_records(roots, graph)
    pending = collect(keys(roots))
    seen = Set{Base.UUID}()
    records = NamedTuple[]
    while !isempty(pending)
        uuid = pop!(pending)
        uuid in seen && continue
        push!(seen, uuid)
        info = get(graph, uuid, nothing)
        info === nothing && return nothing
        record = _dependency_record(uuid, info)
        record === nothing && return nothing
        mod = get(roots, uuid, get(Base.loaded_modules, Base.PkgId(uuid, info.name), nothing))
        mod === nothing ||
            _dependency_module_matches(mod, uuid, info, record) ||
            return nothing
        push!(records, record)
        append!(pending, values(info.dependencies))
    end
    return sort!(records; by = record -> record.uuid)
end

function _julia_dependency_identity(modules; graph = nothing)
    try
        selected = graph === nothing ? Pkg.dependencies() : graph
        roots = Dict(Base.PkgId(mod).uuid => mod for mod in modules)
        records = _dependency_records(roots, selected)
        records === nothing && return nothing
        return (; julia = string(VERSION), packages = records)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _implementation_identity(dependencies = (); source_sha256 = _loaded_source_identity())
    if source_sha256 === nothing
        @warn "Iduna source identity unavailable; affected cache reuse is disabled."
        return nothing
    end
    common = (BioAlignments, BioSequences, MIToS, StatsBase, Statistics, Random,
        CSV, DataFrames, CodecZlib, JSON, HTTP)
    modules = sort!(unique([common..., dependencies...]); by = mod -> String(nameof(mod)))
    dependency_identity = _julia_dependency_identity(modules)
    if dependency_identity === nothing
        @warn "Julia dependency identity unavailable; affected cache reuse is disabled."
        return nothing
    end
    return (; contract = 1, source_sha256, julia = string(VERSION),
        platform = Sys.MACHINE, cpu = Sys.CPU_NAME, cpu_threads = Sys.CPU_THREADS,
        julia_threads = Threads.nthreads(),
        thread_environment = Dict(name => get(ENV, name, nothing)
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "HMMER_NCPU")),
        dependencies = dependency_identity.packages)
end

function _executable_identity(path::AbstractString)
    try
        isfile(path) || return nothing
        resolved = realpath(path)
        return (; path = resolved, sha256 = _file_sha256(resolved))
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _jll_tool_identity(mod::Module, names)
    try
        executables = Dict(String(name) =>
                               _executable_identity(first(getproperty(mod, name)()))
        for name in names)
        any(isnothing, values(executables)) && return nothing
        version = Base.pkgversion(mod)
        version === nothing && return nothing
        isdir(mod.artifact_dir) &&
        occursin(r"^[0-9a-f]{40}$", basename(mod.artifact_dir)) ||
            return nothing
        return (; package = String(nameof(mod)), version = string(version),
            artifact_tree = basename(mod.artifact_dir), executables)
    catch err
        err isa InterruptException && rethrow()
        @warn "External tool identity unavailable; affected cache reuse is disabled." package=String(nameof(mod)) reason=sprint(
            showerror, err)
        return nothing
    end
end
