module IdunaClustalOExt

using ClustalO_jll: ClustalO_jll
using Iduna: Iduna

function Iduna.EPLI.clustalo_aligner(input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        run_label::AbstractString = "run",
        aligner_args::Cmd = Cmd(String[]),
        runner::Function = run,
        aligner = ClustalO_jll.clustalo())
    mkpath(dirname(output_fasta))
    command = `$aligner --infile $input_fasta --outfile $output_fasta --outfmt=fasta --force --output-order=input-order $aligner_args`
    if logs_dir === nothing
        runner(command)
    else
        mkpath(logs_dir)
        stdout_log = joinpath(logs_dir, "$(run_label)_clustalo_stdout.log")
        stderr_log = joinpath(logs_dir, "$(run_label)_clustalo_stderr.log")
        open(stdout_log, "w") do out
            open(stderr_log, "w") do err
                runner(pipeline(command; stdout = out, stderr = err))
            end
        end
    end
    isfile(output_fasta) ||
        error("Clustal Omega did not write the expected MSA FASTA at $(output_fasta).")
    return output_fasta
end

end
