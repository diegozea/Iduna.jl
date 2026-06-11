# Command Execution
# -----------------

function _open_logs(f::Function, stdout_path::AbstractString, stderr_path::AbstractString)
    mkpath(dirname(stdout_path))
    mkpath(dirname(stderr_path))
    open(stdout_path, "w") do out_io
        open(stderr_path, "w") do err_io
            f(out_io, err_io)
        end
    end
end

function _forward_stdout(input::IO, outputs::IO...)
    while !eof(input)
        bytes = readavailable(input)
        isempty(bytes) && continue
        for output in outputs
            write(output, bytes)
            flush(output)
        end
    end
    return nothing
end

function _wait_with_progress(process,
        progress_desc::Union{Nothing, AbstractString};
        progress_output::IO = stderr,
        progress_enabled::Bool = _terminal_progress_enabled(progress_output))
    if progress_desc === nothing || !progress_enabled
        wait(process)
        return nothing
    end

    progress = ProgressMeter.ProgressUnknown(;
        desc = String(progress_desc),
        dt = _TRANSCRIPT_QUERY_SPINNER_INTERVAL_SECONDS,
        spinner = true,
        output = progress_output)
    wait_task = @async wait(process)
    completed = false
    try
        while !istaskdone(wait_task)
            ProgressMeter.next!(progress; force = true)
            timedwait(() -> istaskdone(wait_task),
                _TRANSCRIPT_QUERY_SPINNER_INTERVAL_SECONDS; pollint = 0.1)
        end
        wait(wait_task)
        completed = true
    finally
        if completed && success(process)
            ProgressMeter.finish!(progress)
        else
            ProgressMeter.cancel(progress, "External command stopped")
        end
    end
    return nothing
end

function _run_logged_command(command::Cmd,
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        live_stdout::Bool = false,
        progress_desc::Union{Nothing, AbstractString} = nothing,
        progress_output::IO = stderr,
        progress_enabled::Bool = _terminal_progress_enabled(progress_output))
    _open_logs(stdout_log, stderr_log) do out_io, err_io
        out_pipe = Pipe()
        # When a spinner is shown, keep command output in logs to avoid garbled progress text.
        show_live_stdout = live_stdout && !(progress_desc !== nothing && progress_enabled)
        outputs = show_live_stdout ? (out_io, stderr) : (out_io,)
        stdout_task = @async _forward_stdout(out_pipe, outputs...)
        process = run(pipeline(command; stdout = out_pipe, stderr = err_io), wait = false)
        close(out_pipe.in)
        _wait_with_progress(process, progress_desc;
            progress_output,
            progress_enabled)
        close(out_pipe)
        wait(stdout_task)
        success(process) ||
            error("ThorAxe command failed: $(command). See $(stderr_log).")
    end
    return nothing
end

function _thoraxe_runner(stdout_log::AbstractString,
        stderr_log::AbstractString;
        live_stdout::Bool = false,
        progress_desc::Union{Nothing, AbstractString} = nothing)
    # ThorAxe.jl builds the exact CLI command. The runner only controls logging
    # and error messages around that command.
    return command -> _run_logged_command(
        command, stdout_log, stderr_log; live_stdout, progress_desc)
end

function _transcript_query_thoraxe_runner(stdout_log::AbstractString,
        stderr_log::AbstractString)
    return _thoraxe_runner(stdout_log, stderr_log;
        live_stdout = true,
        progress_desc = "Running ThorAxe transcript_query: ")
end

function _retry_wait_seconds(attempt::Integer)
    return min(30.0, 2.0^(attempt - 1))
end

function _biomart_cache_dir()
    return Scratch.@get_scratch!("biomart_datasets")
end
