"""
    mafft_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run MAFFT on an unaligned FASTA and write the aligned FASTA to `output_fasta`.

This helper is public but intentionally not exported. It is implemented by the
`IdunaMAFFTExt` package extension when `MAFFT_jll` is loaded.
"""
function mafft_aligner end

"""
    clustalo_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run Clustal Omega on an unaligned FASTA and write the aligned FASTA to `output_fasta`.

This helper is public but intentionally not exported. It is implemented by the
`IdunaClustalOExt` package extension when `ClustalO_jll` is loaded.
"""
function clustalo_aligner end
