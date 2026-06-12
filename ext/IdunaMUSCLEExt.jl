module IdunaMUSCLEExt

using Iduna: Iduna
using MUSCLE_jll: MUSCLE_jll

function _read_fasta_records(path::AbstractString)
    records = Tuple{String, String}[]
    current_header = nothing
    current_sequence = String[]
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, '>')
            if current_header !== nothing
                push!(records, (current_header, join(current_sequence)))
            end
            current_header = stripped[2:end]
            empty!(current_sequence)
        else
            current_header === nothing &&
                error("FASTA sequence data appeared before a header in $(path).")
            push!(current_sequence, stripped)
        end
    end
    if current_header !== nothing
        push!(records, (current_header, join(current_sequence)))
    end
    isempty(records) && error("FASTA file has no sequences: $(path).")
    return records
end

function _write_fasta_records(path::AbstractString, records)
    open(path, "w") do io
        for (header, sequence) in records
            println(io, '>', header)
            println(io, sequence)
        end
    end
    return path
end

function _reorder_fasta_like_input!(input_fasta::AbstractString,
        output_fasta::AbstractString)
    input_records = _read_fasta_records(input_fasta)
    output_records = _read_fasta_records(output_fasta)
    length(input_records) == length(output_records) ||
        error("MUSCLE output sequence count does not match the input FASTA.")

    output_by_header = Dict{String, Vector{String}}()
    for (header, sequence) in output_records
        push!(get!(output_by_header, header, String[]), sequence)
    end

    reordered = Tuple{String, String}[]
    for (header, _sequence) in input_records
        sequences = get(output_by_header, header, nothing)
        if sequences === nothing || isempty(sequences)
            error("MUSCLE output is missing FASTA record $(repr(header)).")
        end
        push!(reordered, (header, popfirst!(sequences)))
    end
    if any(sequences -> !isempty(sequences), values(output_by_header))
        error("MUSCLE output contains FASTA records that were not in the input.")
    end

    tmp = tempname(dirname(output_fasta))
    _write_fasta_records(tmp, reordered)
    mv(tmp, output_fasta; force = true)
    return output_fasta
end

function Iduna.EPLI.muscle_aligner(input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        run_label::AbstractString = "run",
        aligner_args::Cmd = Cmd(String[]),
        runner::Function = run,
        aligner = MUSCLE_jll.muscle())
    mkpath(dirname(output_fasta))
    command = `$aligner -align $input_fasta -output $output_fasta $aligner_args`
    if logs_dir === nothing
        runner(command)
    else
        mkpath(logs_dir)
        stdout_log = joinpath(logs_dir, "$(run_label)_muscle_stdout.log")
        stderr_log = joinpath(logs_dir, "$(run_label)_muscle_stderr.log")
        open(stdout_log, "w") do out
            open(stderr_log, "w") do err
                runner(pipeline(command; stdout = out, stderr = err))
            end
        end
    end
    isfile(output_fasta) ||
        error("MUSCLE did not write the expected MSA FASTA at $(output_fasta).")
    _reorder_fasta_like_input!(input_fasta, output_fasta)
    return output_fasta
end

end
