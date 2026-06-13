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
    famsa_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run FAMSA on an unaligned FASTA and write the aligned FASTA to `output_fasta`.
Iduna rewrites the aligned FASTA records in the same order as `input_fasta`.
This row-order rewrite uses MIToS FASTA parsing, so sequences should use the
MIToS residue alphabet.

# Fixed Arguments

Iduna runs FAMSA with user `aligner_args`, followed by the positional
`input_fasta output_fasta` arguments.

This helper is public but intentionally not exported. It is implemented by the
`IdunaFAMSAExt` package extension when `FAMSA_jll` is loaded.
"""
function famsa_aligner end

"""
    kalign_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run Kalign on an unaligned FASTA and write the aligned FASTA to `output_fasta`.
Iduna rewrites the aligned FASTA records in the same order as `input_fasta`.
This row-order rewrite uses MIToS FASTA parsing, so sequences should use the
MIToS residue alphabet.

# Fixed Arguments

Iduna runs Kalign with `-i input_fasta -o output_fasta --format fasta`, then
appends `aligner_args`. Use `aligner_args` to pass `--type` when you want
Kalign to validate a specific sequence type.

This helper is public but intentionally not exported. It is implemented by the
`IdunaKalignExt` package extension when `kalign_jll` is loaded.
"""
function kalign_aligner end

"""
    muscle_aligner(input_fasta, output_fasta; logs_dir=nothing, run_label="run",
        aligner_args=Cmd(String[]))

Run MUSCLE on an unaligned FASTA and write the aligned FASTA to `output_fasta`.
Iduna rewrites the aligned FASTA records in the same order as `input_fasta`,
because MUSCLE may emit guide-tree order by default.
This row-order rewrite uses MIToS FASTA parsing, so sequences should use the
MIToS residue alphabet.

# Fixed Arguments

Iduna runs MUSCLE with `-align input_fasta -output output_fasta`, then appends
`aligner_args`.

This helper is public but intentionally not exported. It is implemented by the
`IdunaMUSCLEExt` package extension when `MUSCLE_jll` is loaded.
"""
function muscle_aligner end

function _read_fasta_records(path::AbstractString)
    records = [(String(sequence_id(sequence)), sequence)
               for sequence in read_file(path, FASTASequences)]
    isempty(records) && error("FASTA file has no sequences: $(path).")
    return records
end

function _write_fasta_records(path::AbstractString,
        msa::AbstractMultipleSequenceAlignment)
    write_file(path, msa, FASTA)
    return path
end

function _write_fasta_records(path::AbstractString, records)
    sequences = [sequence isa AbstractSequence ?
                 AnnotatedSequence(sequence) :
                 AnnotatedSequence(String(header), String(sequence))
                 for (header, sequence) in records]
    write_file(path, sequences, FASTASequences)
    return path
end

function _reorder_fasta_like_input!(input_fasta::AbstractString,
        output_fasta::AbstractString;
        aligner_name::AbstractString = "Aligner")
    input_records = _read_fasta_records(input_fasta)
    input_names = [name for (name, _sequence) in input_records]
    output_msa = read_file(output_fasta, FASTA)
    nsequences(output_msa) == length(input_names) ||
        error("$(aligner_name) output sequence count does not match the input FASTA.")

    output_names = Set(String.(sequencenames(output_msa)))
    for name in input_names
        if !(name in output_names)
            error("$(aligner_name) output is missing FASTA record $(repr(name)).")
        end
    end
    if any(name -> !(name in input_names), String.(sequencenames(output_msa)))
        error("$(aligner_name) output contains FASTA records that were not in the input.")
    end

    reordered = output_msa[input_names, :]
    output_dir = dirname(output_fasta)
    tmp = tempname(isempty(output_dir) ? "." : output_dir)
    # This intentionally writes through MIToS FASTA; bundled JLL aligners are
    # documented for MIToS-compatible residue alphabets.
    _write_fasta_records(tmp, reordered)
    mv(tmp, output_fasta; force = true)
    return output_fasta
end
