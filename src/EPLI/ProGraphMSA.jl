"""
    prographmsa_aligner(input_fasta, output_fasta; logs_dir=nothing, sample_label="sample")

Run ProGraphMSA on an unaligned FASTA and write the aligned FASTA to `output_fasta`.

This helper is public but intentionally not exported because it wraps separate
software. Use it as `Iduna.EPLI.prographmsa_aligner`.
"""
function prographmsa_aligner(input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        sample_label::AbstractString = "sample",
        runner::Function = run,
        aligner::Union{Nothing, AbstractString} = nothing)
    aligner_path = aligner === nothing ? ThorAxe.aligner_executable() : String(aligner)
    mkpath(dirname(output_fasta))
    command = `$aligner_path --input_order --fasta --output $output_fasta $input_fasta`
    if logs_dir === nothing
        runner(command)
    else
        mkpath(logs_dir)
        stdout_log = joinpath(logs_dir, "$(sample_label)_prographmsa_stdout.log")
        stderr_log = joinpath(logs_dir, "$(sample_label)_prographmsa_stderr.log")
        open(stdout_log, "w") do out
            open(stderr_log, "w") do err
                runner(pipeline(command; stdout = out, stderr = err))
            end
        end
    end
    isfile(output_fasta) ||
        error("ProGraphMSA did not write the expected MSA FASTA at $(output_fasta).")
    return output_fasta
end
