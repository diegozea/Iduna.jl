"""
    mafft_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run MAFFT on an unaligned FASTA and write the aligned FASTA to `output_fasta`.
Iduna always keeps aligned FASTA records in the same order as `input_fasta`.

# Fixed Arguments

Iduna runs MAFFT with `--quiet --inputorder`, then appends `aligner_args` before
the input FASTA path.

This helper is public but intentionally not exported. It is implemented by the
`IdunaMAFFTExt` package extension when `MAFFT_jll` is loaded.
"""
function mafft_aligner end

"""
    clustalo_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run Clustal Omega on an unaligned FASTA and write the aligned FASTA to
`output_fasta`. Iduna always keeps aligned FASTA records in the same order as
`input_fasta`.

# Fixed Arguments

Iduna runs Clustal Omega with `--infile input_fasta --outfile output_fasta
--outfmt=fasta --force --output-order=input-order`, then appends
`aligner_args`.

This helper is public but intentionally not exported. It is implemented by the
`IdunaClustalOExt` package extension when `ClustalO_jll` is loaded.
"""
function clustalo_aligner end

"""
    muscle_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run MUSCLE on an unaligned FASTA and write the aligned FASTA to `output_fasta`.
Iduna rewrites the aligned FASTA records in the same order as `input_fasta`,
because MUSCLE may emit guide-tree order by default.

# Fixed Arguments

Iduna runs MUSCLE with `-align input_fasta -output output_fasta`, then appends
`aligner_args`.

This helper is public but intentionally not exported. It is implemented by the
`IdunaMUSCLEExt` package extension when `MUSCLE_jll` is loaded.
"""
function muscle_aligner end
