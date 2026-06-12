"""
    prographmsa_aligner(input_fasta, output_fasta; logs_dir=nothing,
        run_label="run", aligner_args=Cmd(String[]))

Run ProGraphMSA on an unaligned FASTA and write the aligned FASTA to
`output_fasta`. Iduna always keeps aligned FASTA records in the same order as
`input_fasta`.

# Fixed Arguments

Iduna runs ProGraphMSA with `--input_order --fasta --output output_fasta`, then
appends `aligner_args` before the input FASTA path.

This helper is public but intentionally not exported because it wraps separate
software. Use it as `Iduna.EPLI.prographmsa_aligner`.
"""
function prographmsa_aligner(input_fasta::AbstractString,
        output_fasta::AbstractString;
        logs_dir::Union{Nothing, AbstractString} = nothing,
        run_label::AbstractString = "run",
        aligner_args::Cmd = Cmd(String[]),
        runner::Function = run,
        aligner::Union{Nothing, AbstractString} = nothing)
    aligner_path = aligner === nothing ? ThorAxe.aligner_executable() : String(aligner)
    mkpath(dirname(output_fasta))
    command = `$aligner_path --input_order --fasta --output $output_fasta $aligner_args $input_fasta`
    if logs_dir === nothing
        runner(command)
    else
        mkpath(logs_dir)
        stdout_log = joinpath(logs_dir, "$(run_label)_prographmsa_stdout.log")
        stderr_log = joinpath(logs_dir, "$(run_label)_prographmsa_stderr.log")
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
