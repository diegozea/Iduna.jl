module IdunaFAMSAExt

using FAMSA_jll: FAMSA_jll
using Iduna: Iduna

function Iduna.EPLI.famsa_aligner(input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        run_label::AbstractString = "run",
        aligner_args::Cmd = Cmd(String[]),
        runner::Function = run,
        aligner = FAMSA_jll.famsa())
    mkpath(dirname(output_fasta))
    command = `$aligner $aligner_args $input_fasta $output_fasta`
    if logs_dir === nothing
        runner(command)
    else
        mkpath(logs_dir)
        stdout_log = joinpath(logs_dir, "$(run_label)_famsa_stdout.log")
        stderr_log = joinpath(logs_dir, "$(run_label)_famsa_stderr.log")
        open(stdout_log, "w") do out
            open(stderr_log, "w") do err
                runner(pipeline(command; stdout = out, stderr = err))
            end
        end
    end
    isfile(output_fasta) ||
        error("FAMSA did not write the expected MSA FASTA at $(output_fasta).")
    Iduna.EPLI._reorder_fasta_like_input!(input_fasta, output_fasta;
        aligner_name = "FAMSA")
    return output_fasta
end

end
