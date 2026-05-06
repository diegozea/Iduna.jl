module ResultsValidation

import CSV

using DataFrames: DataFrame
using MIToS.MSA: A3M, AbstractMultipleSequenceAlignment, FASTA, Stockholm,
                 hobohmI, n_effective, namedmatrix, ncolumns, nsequences,
                 read_file, sequencenames, stringsequence

using ..Utils: ExpansionResult, ResolvedTarget, SeedSelection, ValidationResult,
               protein_alignment_stats, resolve_sequence_name

export alignment_stats,
       load_expanded_msa,
       load_seed_msa,
       load_msa,
       validate_results

function load_msa(path::AbstractString; keepinserts::Bool = true)
    lower = lowercase(path)
    if endswith(lower, ".sto") || endswith(lower, ".stockholm")
        return read_file(path, Stockholm; keepinserts)
    elseif endswith(lower, ".a3m")
        return read_file(path, A3M; keepinserts)
    elseif endswith(lower, ".fasta") || endswith(lower, ".fa")
        return read_file(path, FASTA)
    end
    error("Cannot infer MSA format from $(path).")
end

function load_seed_msa(seed::SeedSelection; keepinserts::Bool = true)
    load_msa(seed.stockholm_path; keepinserts)
end

function load_expanded_msa(expansion::ExpansionResult; keepinserts::Bool = true)
    load_msa(expansion.match_stockholm; keepinserts)
end

function alignment_stats(path::AbstractString; cluster_threshold::Real = 62.0, neff_threshold::Real = 80.0)
    msa = load_msa(path; keepinserts = true)
    # Hobohm clustering gives a compact diversity summary for the alignment.
    clusters = hobohmI(namedmatrix(msa), cluster_threshold)
    return (;
        msa,
        n_sequences = nsequences(msa),
        n_columns = ncolumns(msa),
        n_clusters = length(clusters.clustersize),
        neff = Float64(n_effective(msa, neff_threshold))
    )
end

function _resolve_query_name(msa::AbstractMultipleSequenceAlignment,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    return resolve_sequence_name(msa, (gene_id, transcript_id); fallback = true)
end

function _extract_query_sequence(msa::AbstractMultipleSequenceAlignment,
        gene_id::AbstractString,
        transcript_id::AbstractString)
    name = _resolve_query_name(msa, gene_id, transcript_id)
    seq = replace(stringsequence(msa, name), '-' => "", '.' => "")
    return name, uppercase(String(seq))
end

function _read_fasta_sequence(path::AbstractString)
    msa = read_file(path, FASTA)
    name = first(sequencenames(msa))
    return uppercase(String(stringsequence(msa, name)))
end

function _align_sequences(query_seq::AbstractString, reference_seq::AbstractString)
    protein_alignment_stats(query_seq, reference_seq; include_alignment = true)
end

function _write_alignment_log(path::AbstractString,
        target::ResolvedTarget,
        query_name::AbstractString,
        query_seq::AbstractString,
        reference_seq::AbstractString,
        aln_stats)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Query vs UniProt alignment")
        println(io, "input_id=$(target.input_id)")
        println(io, "uniprot_id=$(target.uniprot_id)")
        println(io, "gene_id=$(target.ensembl_gene_id)")
        println(io, "transcript_id=$(target.transcript_id)")
        println(io, "query_name=$(query_name)")
        println(io, "query_length=$(length(query_seq))")
        println(io, "uniprot_length=$(length(reference_seq))")
        println(io, "mismatches=$(aln_stats.mismatches)")
        println(io, "insertions=$(aln_stats.insertions)")
        println(io, "deletions=$(aln_stats.deletions)")
        println(io)
        show(io, MIME("text/plain"), aln_stats.aln)
        println(io)
    end
    return path
end

function validate_results(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::ExpansionResult,
        workdir::AbstractString)
    seed_stats = alignment_stats(seed.stockholm_path)
    expanded_stats = alignment_stats(expansion.match_stockholm)

    query_name = nothing
    aln_path = nothing
    aln_identical = nothing
    aln_mismatches = nothing
    aln_insertions = nothing
    aln_deletions = nothing
    warnings = String[]

    if target.uniprot_sequence_path !== nothing && isfile(target.uniprot_sequence_path)
        # Compare the final query sequence against UniProt when a reference exists.
        query_name,
        query_seq = _extract_query_sequence(expanded_stats.msa,
            target.ensembl_gene_id, target.transcript_id)
        uniprot_seq = _read_fasta_sequence(target.uniprot_sequence_path)
        aln_stats = _align_sequences(query_seq, uniprot_seq)
        aln_path = joinpath(workdir, "validation", "query_vs_uniprot_alignment.txt")
        _write_alignment_log(
            aln_path, target, query_name, query_seq, uniprot_seq, aln_stats)
        aln_identical = aln_stats.identical
        aln_mismatches = aln_stats.mismatches
        aln_insertions = aln_stats.insertions
        aln_deletions = aln_stats.deletions
        if aln_insertions != 0 || aln_deletions != 0
            push!(warnings, "Expanded query has indels relative to the UniProt sequence.")
        elseif aln_identical == false
            push!(warnings, "Expanded query has substitutions relative to the UniProt sequence.")
        end
    end

    stats_path = joinpath(workdir, "validation", "stats.csv")
    mkpath(dirname(stats_path))
    CSV.write(stats_path,
        DataFrame(
            input_id = [target.input_id],
            uniprot_id = [target.uniprot_id],
            gene_id = [target.ensembl_gene_id],
            transcript_id = [target.transcript_id],
            pid = [seed.pid],
            seed_path = [seed.stockholm_path],
            expanded_path = [expansion.match_stockholm],
            seed_nseq = [seed_stats.n_sequences],
            expanded_nseq = [expanded_stats.n_sequences],
            seed_ncol = [seed_stats.n_columns],
            expanded_ncol = [expanded_stats.n_columns],
            seed_clusters62 = [seed_stats.n_clusters],
            expanded_clusters62 = [expanded_stats.n_clusters],
            seed_neff80 = [seed_stats.neff],
            expanded_neff80 = [expanded_stats.neff],
            query_name = [query_name],
            query_vs_uniprot_path = [aln_path],
            aln_identical = [aln_identical],
            aln_mismatches = [aln_mismatches],
            aln_insertions = [aln_insertions],
            aln_deletions = [aln_deletions]
        );
        transform = (_, val) -> something(val, missing))

    status = isempty(warnings) ? :ok : :warn
    return ValidationResult(;
        stats_path,
        query_name,
        query_vs_uniprot_path = aln_path,
        seed_nseq = seed_stats.n_sequences,
        seed_ncol = seed_stats.n_columns,
        seed_clusters62 = seed_stats.n_clusters,
        seed_neff80 = seed_stats.neff,
        expanded_nseq = expanded_stats.n_sequences,
        expanded_ncol = expanded_stats.n_columns,
        expanded_clusters62 = expanded_stats.n_clusters,
        expanded_neff80 = expanded_stats.neff,
        aln_identical,
        aln_mismatches,
        aln_insertions,
        aln_deletions,
        warnings,
        status
    )
end

end
