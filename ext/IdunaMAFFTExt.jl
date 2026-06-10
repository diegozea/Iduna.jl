module IdunaMAFFTExt

using Iduna: Iduna
using MAFFT_jll: MAFFT_jll

function Iduna.EPLI.mafft_aligner(input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        run_label::AbstractString = "run",
        aligner_args::Cmd = Cmd(String[]),
        runner::Function = run,
        aligner = MAFFT_jll.mafft())
    mkpath(dirname(output_fasta))
    command = `$aligner --quiet --inputorder $aligner_args $input_fasta`
    if logs_dir === nothing
        runner(pipeline(command; stdout = output_fasta))
    else
        mkpath(logs_dir)
        stderr_log = joinpath(logs_dir, "$(run_label)_mafft_stderr.log")
        open(stderr_log, "w") do err
            runner(pipeline(command; stdout = output_fasta, stderr = err))
        end
    end
    isfile(output_fasta) ||
        error("MAFFT did not write the expected MSA FASTA at $(output_fasta).")
    return output_fasta
end

end
