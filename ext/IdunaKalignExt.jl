module IdunaKalignExt

using Iduna: Iduna
using kalign_jll: kalign_jll

function Iduna.EPLI.kalign_aligner(input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        run_label::AbstractString = "run",
        aligner_args::Cmd = Cmd(String[]),
        runner::Function = run,
        aligner = kalign_jll.kalign())
    mkpath(dirname(output_fasta))
    # Kalign validates --type against its own autodetected alphabet. Do not force
    # protein mode here: proteins containing only A/C/G/T are valid EPLI inputs,
    # but Kalign detects them as DNA and rejects --type protein before alignment.
    command = `$aligner -i $input_fasta -o $output_fasta --format fasta $aligner_args`
    if logs_dir === nothing
        runner(command)
    else
        mkpath(logs_dir)
        stdout_log = joinpath(logs_dir, "$(run_label)_kalign_stdout.log")
        stderr_log = joinpath(logs_dir, "$(run_label)_kalign_stderr.log")
        open(stdout_log, "w") do out
            open(stderr_log, "w") do err
                runner(pipeline(command; stdout = out, stderr = err))
            end
        end
    end
    isfile(output_fasta) ||
        error("Kalign did not write the expected MSA FASTA at $(output_fasta).")
    Iduna.EPLI._reorder_fasta_like_input!(input_fasta, output_fasta;
        aligner_name = "Kalign")
    return output_fasta
end

end
