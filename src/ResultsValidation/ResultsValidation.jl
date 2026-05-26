module ResultsValidation

import CSV

using DataFrames: DataFrame
using MIToS.MSA: A3M, AbstractMultipleSequenceAlignment, FASTA, Stockholm,
                 hobohmI, n_effective, namedmatrix, ncolumns, nsequences,
                 read_file, sequencenames, stringsequence

using ..Utils: ExpansionResult, ResolvedTarget, SeedSelection, ValidationResult,
               _relative_artifact_path, _resolve_artifact_path, format_pid_dir,
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

function load_seed_msa(seed::SeedSelection; keepinserts::Bool = true,
        workdir::Union{Nothing, AbstractString} = seed.workdir)
    path = workdir === nothing ? seed.stockholm_path :
           _resolve_artifact_path(seed.stockholm_path, workdir)
    load_msa(path; keepinserts)
end

function load_expanded_msa(expansion::ExpansionResult; keepinserts::Bool = true,
        workdir::Union{Nothing, AbstractString} = expansion.workdir)
    path = workdir === nothing ? expansion.match_stockholm :
           _resolve_artifact_path(expansion.match_stockholm, workdir)
    load_msa(path; keepinserts)
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

function _validation_input_paths(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::Union{Nothing, ExpansionResult},
        workdir::AbstractString)
    seed_path = _resolve_artifact_path(
        seed.stockholm_path, seed.workdir === nothing ? workdir : seed.workdir)
    expanded_path = expansion === nothing ? nothing :
                    _resolve_artifact_path(expansion.match_stockholm,
        expansion.workdir === nothing ? workdir : expansion.workdir)
    uniprot_path = target.uniprot_sequence_path === nothing ? nothing :
                   _resolve_artifact_path(target.uniprot_sequence_path,
        target.workdir === nothing ? workdir : target.workdir)
    return (; seed_path, expanded_path, uniprot_path)
end

function _validation_alignment_stats(paths, expansion::Union{Nothing, ExpansionResult})
    seed_stats = alignment_stats(paths.seed_path)
    expanded_stats = expansion === nothing ? nothing :
                     alignment_stats(paths.expanded_path)
    final_stats = expansion === nothing ? seed_stats : expanded_stats
    return (; seed_stats, expanded_stats, final_stats)
end

function _empty_uniprot_comparison()
    return (;
        query_name = nothing,
        aln_path = nothing,
        aln_identical = nothing,
        aln_mismatches = nothing,
        aln_insertions = nothing,
        aln_deletions = nothing,
        warnings = String[]
    )
end

function _alignment_warnings(expansion::Union{Nothing, ExpansionResult}, aln_stats)
    label = expansion === nothing ? "Seed query" : "Expanded query"
    if aln_stats.insertions != 0 || aln_stats.deletions != 0
        return ["$(label) has indels relative to the UniProt sequence."]
    elseif aln_stats.identical == false
        return ["$(label) has substitutions relative to the UniProt sequence."]
    end
    return String[]
end

function _compare_final_query_to_uniprot(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::Union{Nothing, ExpansionResult},
        final_stats,
        uniprot_path::Union{Nothing, AbstractString},
        workdir::AbstractString)
    (uniprot_path !== nothing && isfile(uniprot_path)) ||
        return _empty_uniprot_comparison()

    # Compare the final available query sequence against UniProt when a reference exists.
    query_name,
    query_seq = _extract_query_sequence(
        final_stats.msa, target.ensembl_gene_id, target.transcript_id)
    uniprot_seq = _read_fasta_sequence(uniprot_path)
    aln_stats = _align_sequences(query_seq, uniprot_seq)
    validation_dir = joinpath(workdir, "validation", format_pid_dir(seed.pid))
    aln_path = joinpath(validation_dir, "query_vs_uniprot_alignment.txt")
    _write_alignment_log(aln_path, target, query_name, query_seq, uniprot_seq, aln_stats)
    return (;
        query_name,
        aln_path,
        aln_identical = aln_stats.identical,
        aln_mismatches = aln_stats.mismatches,
        aln_insertions = aln_stats.insertions,
        aln_deletions = aln_stats.deletions,
        warnings = _alignment_warnings(expansion, aln_stats)
    )
end

function _write_validation_stats(target::ResolvedTarget,
        seed::SeedSelection,
        stats,
        comparison,
        paths,
        workdir::AbstractString)
    stats_path = joinpath(workdir, "validation", format_pid_dir(seed.pid), "stats.csv")
    mkpath(dirname(stats_path))
    seed_summary_path = _relative_artifact_path(paths.seed_path, workdir)
    expanded_summary_path = _relative_artifact_path(paths.expanded_path, workdir)
    aln_summary_path = _relative_artifact_path(comparison.aln_path, workdir)
    seed_stats = stats.seed_stats
    expanded_stats = stats.expanded_stats
    expanded_nseq = expanded_stats === nothing ? nothing : expanded_stats.n_sequences
    expanded_ncol = expanded_stats === nothing ? nothing : expanded_stats.n_columns
    expanded_clusters62 = expanded_stats === nothing ? nothing : expanded_stats.n_clusters
    expanded_neff80 = expanded_stats === nothing ? nothing : expanded_stats.neff
    CSV.write(stats_path,
        DataFrame(
            input_id = [target.input_id],
            uniprot_id = [target.uniprot_id],
            gene_id = [target.ensembl_gene_id],
            transcript_id = [target.transcript_id],
            pid = [seed.pid],
            seed_path = [seed_summary_path],
            expanded_path = [expanded_summary_path],
            seed_nseq = [seed_stats.n_sequences],
            expanded_nseq = [expanded_nseq],
            seed_ncol = [seed_stats.n_columns],
            expanded_ncol = [expanded_ncol],
            seed_clusters62 = [seed_stats.n_clusters],
            expanded_clusters62 = [expanded_clusters62],
            seed_neff80 = [seed_stats.neff],
            expanded_neff80 = [expanded_neff80],
            query_name = [comparison.query_name],
            query_vs_uniprot_path = [aln_summary_path],
            aln_identical = [comparison.aln_identical],
            aln_mismatches = [comparison.aln_mismatches],
            aln_insertions = [comparison.aln_insertions],
            aln_deletions = [comparison.aln_deletions]
        );
        transform = (_, val) -> something(val, missing))
    return stats_path
end

function _validation_result(stats_path::AbstractString, stats, comparison, warnings)
    seed_stats = stats.seed_stats
    expanded_stats = stats.expanded_stats
    expanded_nseq = expanded_stats === nothing ? nothing : expanded_stats.n_sequences
    expanded_ncol = expanded_stats === nothing ? nothing : expanded_stats.n_columns
    expanded_clusters62 = expanded_stats === nothing ? nothing : expanded_stats.n_clusters
    expanded_neff80 = expanded_stats === nothing ? nothing : expanded_stats.neff
    status = isempty(warnings) ? :ok : :warn
    return ValidationResult(;
        stats_path,
        query_name = comparison.query_name,
        query_vs_uniprot_path = comparison.aln_path,
        seed_nseq = seed_stats.n_sequences,
        seed_ncol = seed_stats.n_columns,
        seed_clusters62 = seed_stats.n_clusters,
        seed_neff80 = seed_stats.neff,
        expanded_nseq,
        expanded_ncol,
        expanded_clusters62,
        expanded_neff80,
        aln_identical = comparison.aln_identical,
        aln_mismatches = comparison.aln_mismatches,
        aln_insertions = comparison.aln_insertions,
        aln_deletions = comparison.aln_deletions,
        warnings,
        status
    )
end

function validate_results(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::Union{Nothing, ExpansionResult},
        workdir::AbstractString)
    paths = _validation_input_paths(target, seed, expansion, workdir)
    stats = _validation_alignment_stats(paths, expansion)
    comparison = _compare_final_query_to_uniprot(
        target, seed, expansion, stats.final_stats, paths.uniprot_path, workdir)
    warnings = comparison.warnings
    stats_path = _write_validation_stats(target, seed, stats, comparison, paths, workdir)
    return _validation_result(stats_path, stats, comparison, warnings)
end

end
