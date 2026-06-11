"""
    format_pid(pid) -> String

Format a percent identity value for user-facing tables and messages.

# Arguments

- `pid::Real`: percent identity value.
"""
format_pid(pid::Real) = @sprintf("%.1f", Float64(pid))

"""
    format_pid_dir(pid) -> String

Format a percent identity value as a stable directory name.

# Arguments

- `pid::Real`: percent identity value.
"""
format_pid_dir(pid::Real) = "pid_$(@sprintf("%.2f", Float64(pid)))"

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

"""
    safe_rm(path, root)

Remove a path only when it is inside the active work directory.

# Arguments

- `path::AbstractString`: file or directory to remove.
- `root::AbstractString`: work directory that must contain `path`.

# Throws

- `ErrorException`: if `path` is outside `root`.
"""
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
