"""
    ResultsValidation

Load MSAs and summarize seed and expanded alignment quality for Iduna results.
"""
module ResultsValidation

import CSV

using DataFrames: DataFrame
using MIToS.MSA: A3M, AbstractMultipleSequenceAlignment, FASTA, Stockholm,
                 hobohmI, n_effective, namedmatrix, ncolumns, nsequences,
                 read_file, sequencenames, stringsequence

using ..Utils: ExpansionResult, ResolvedTarget, SeedSelection, ValidationResult,
               _relative_artifact_path, _resolve_artifact_path, format_pid_dir,
               protein_alignment_stats, resolve_sequence_name, safe_rm,
               _classify_stage_state, _file_sha256, _read_stage_state,
               _write_stage_state

# Expansion Presence
# ------------------

const _MaybeExpansion = Union{Nothing, Missing, ExpansionResult}

function _has_expansion(expansion::_MaybeExpansion)
    !(expansion === nothing || ismissing(expansion))
end

export alignment_stats,
       load_expanded_msa,
       load_seed_msa,
       load_msa,
       validate_results

# MSA Loading
# -----------

"""
    load_msa(path; keepinserts=true)

Load an MSA from Stockholm, A3M, FASTA, or `.fa` input.

# Arguments

- `path::AbstractString`: MSA file path to load.

# Keywords

- `keepinserts::Bool = true`: keep insertion columns when reading formats that
  support this choice.
"""
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

"""
    load_seed_msa(seed; keepinserts=true, workdir=seed.workdir)

Load the selected ThorAxe seed MSA described by a [`SeedSelection`](@ref).

# Arguments

- `seed::SeedSelection`: selected ThorAxe seed to load.

# Keywords

- `keepinserts::Bool = true`: keep insertion columns when reading the MSA.
- `workdir = seed.workdir`: work directory used to resolve relative result
  paths.
"""
function load_seed_msa(seed::SeedSelection; keepinserts::Bool = true,
        workdir::Union{Nothing, AbstractString} = seed.workdir)
    path = workdir === nothing ? seed.stockholm_path :
           _resolve_artifact_path(seed.stockholm_path, workdir)
    load_msa(path; keepinserts)
end

"""
    load_expanded_msa(expansion; keepinserts=true, workdir=expansion.workdir)

Load the match-column expanded MSA described by an [`ExpansionResult`](@ref).

# Arguments

- `expansion::ExpansionResult`: expansion result to load.

# Keywords

- `keepinserts::Bool = true`: keep insertion columns when reading the MSA.
- `workdir = expansion.workdir`: work directory used to resolve relative result
  paths.
"""
function load_expanded_msa(expansion::ExpansionResult; keepinserts::Bool = true,
        workdir::Union{Nothing, AbstractString} = expansion.workdir)
    path = workdir === nothing ? expansion.match_stockholm :
           _resolve_artifact_path(expansion.match_stockholm, workdir)
    load_msa(path; keepinserts)
end

"""
    alignment_stats(path; cluster_threshold=62.0, neff_threshold=80.0)

Load one MSA and compute simple size and diversity statistics.

# Arguments

- `path::AbstractString`: MSA file path to load and summarize.

# Keywords

- `cluster_threshold::Real = 62.0`: identity threshold for Hobohm clustering.
- `neff_threshold::Real = 80.0`: identity threshold for effective sequence count.
  Both thresholds are percent identity values.

# Returns

- A named tuple with the MSA, number of sequences, number of columns, number of
  clusters, and effective sequence count.
"""
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

# Query Sequence Comparison
# -------------------------

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

# Validation State
# ----------------

function _validation_input_paths(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::_MaybeExpansion,
        workdir::AbstractString)
    seed_path = _resolve_artifact_path(
        seed.stockholm_path, seed.workdir === nothing ? workdir : seed.workdir)
    expanded_path = _has_expansion(expansion) ?
                    _resolve_artifact_path(expansion.match_stockholm,
        expansion.workdir === nothing ? workdir : expansion.workdir) : nothing
    uniprot_path = target.uniprot_sequence_path === nothing ? nothing :
                   _resolve_artifact_path(target.uniprot_sequence_path,
        target.workdir === nothing ? workdir : target.workdir)
    return (; seed_path, expanded_path, uniprot_path)
end

function _validation_dir(workdir::AbstractString, seed::SeedSelection)
    joinpath(workdir,
        "validation", format_pid_dir(seed.pid))
end

function _validation_outputs(workdir::AbstractString,
        seed::SeedSelection,
        paths)
    validation_dir = _validation_dir(workdir, seed)
    return (;
        stats_path = joinpath(validation_dir, "stats.csv"),
        query_vs_uniprot_path = paths.uniprot_path === nothing ||
                                !isfile(paths.uniprot_path) ? nothing :
                                joinpath(validation_dir,
            "query_vs_uniprot_alignment.txt"))
end

function _hash_existing(path::Union{Nothing, AbstractString})
    path === nothing && return nothing
    isfile(path) || return nothing
    return _file_sha256(path)
end

function _validation_identity(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::_MaybeExpansion,
        paths)
    return (;
        target = (;
            input_id = target.input_id,
            input_kind = String(target.input_kind),
            uniprot_id = target.uniprot_id,
            ensembl_gene_id = target.ensembl_gene_id,
            transcript_id = target.transcript_id,
            ensembl_protein_id = target.ensembl_protein_id,
            species = target.species),
        seed = (;
            pid = Float64(seed.pid),
            stockholm_sha256 = _hash_existing(paths.seed_path)),
        expansion = _has_expansion(expansion) ?
                    (; match_stockholm_sha256 = _hash_existing(paths.expanded_path),) :
                    nothing,
        # File hashes catch changed alignments even when the path stays the same.
        uniprot_sequence_sha256 = _hash_existing(paths.uniprot_path)
    )
end

function _write_validation_state(workdir::AbstractString,
        seed::SeedSelection,
        status::Symbol,
        identity,
        outputs;
        action,
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing)
    return _write_stage_state(_validation_dir(workdir, seed);
        stage = "validation",
        stage_key = "validation:$(format_pid_dir(seed.pid))",
        status,
        identity,
        outputs,
        warnings,
        exception,
        action,
        workdir)
end

_nothing_if_missing(value) = value === missing ? nothing : value
_maybe_int(value) = value === missing ? nothing : Int(value)
_maybe_float(value) = value === missing ? nothing : Float64(value)
_maybe_bool(value) = value === missing ? nothing : Bool(value)

function _cached_validation_warnings(workdir::AbstractString, seed::SeedSelection)
    state = _read_stage_state(_validation_dir(workdir, seed))
    state isa AbstractDict || return String[]
    return String.(get(state, "warnings", String[]))
end

function _cached_alignment_warnings(row, expansion::_MaybeExpansion)
    identical = _maybe_bool(row.aln_identical)
    identical === nothing && return String[]
    insertions = something(_maybe_int(row.aln_insertions), 0)
    deletions = something(_maybe_int(row.aln_deletions), 0)
    label = _has_expansion(expansion) ? "Expanded query" : "Seed query"
    if insertions != 0 || deletions != 0
        return ["$(label) has indels relative to the UniProt sequence."]
    elseif identical == false
        return ["$(label) has substitutions relative to the UniProt sequence."]
    end
    return String[]
end

function _cached_validation_result(outputs, workdir::AbstractString, seed::SeedSelection,
        expansion::_MaybeExpansion)
    df = DataFrame(CSV.File(outputs.stats_path))
    isempty(df) && error("Cached validation stats at $(outputs.stats_path) are empty.")
    row = first(eachrow(df))
    # Recreate warnings from saved stats so cached and fresh results report the same status.
    warnings = unique(vcat(
        _cached_validation_warnings(workdir, seed),
        _cached_alignment_warnings(row, expansion)))
    return ValidationResult(;
        stats_path = outputs.stats_path,
        query_name = _nothing_if_missing(row.query_name),
        query_vs_uniprot_path = outputs.query_vs_uniprot_path,
        seed_nseq = _maybe_int(row.seed_nseq),
        seed_ncol = _maybe_int(row.seed_ncol),
        seed_clusters62 = _maybe_int(row.seed_clusters62),
        seed_neff80 = _maybe_float(row.seed_neff80),
        expanded_nseq = _maybe_int(row.expanded_nseq),
        expanded_ncol = _maybe_int(row.expanded_ncol),
        expanded_clusters62 = _maybe_int(row.expanded_clusters62),
        expanded_neff80 = _maybe_float(row.expanded_neff80),
        aln_identical = _maybe_bool(row.aln_identical),
        aln_mismatches = _maybe_int(row.aln_mismatches),
        aln_insertions = _maybe_int(row.aln_insertions),
        aln_deletions = _maybe_int(row.aln_deletions),
        warnings,
        status = isempty(warnings) ? :ok : :warn)
end

# Validation Statistics
# ---------------------

function _validation_alignment_stats(paths, expansion::_MaybeExpansion)
    seed_stats = alignment_stats(paths.seed_path)
    expanded_stats = _has_expansion(expansion) ? alignment_stats(paths.expanded_path) :
                     nothing
    final_stats = _has_expansion(expansion) ? expanded_stats : seed_stats
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

function _alignment_warnings(expansion::_MaybeExpansion, aln_stats)
    label = _has_expansion(expansion) ? "Expanded query" : "Seed query"
    if aln_stats.insertions != 0 || aln_stats.deletions != 0
        return ["$(label) has indels relative to the UniProt sequence."]
    elseif aln_stats.identical == false
        return ["$(label) has substitutions relative to the UniProt sequence."]
    end
    return String[]
end

function _compare_final_query_to_uniprot(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::_MaybeExpansion,
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

# Result Construction
# -------------------

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

# Public API
# ----------

"""
    validate_results(target, seed, expansion, workdir; overwrite=false)

Compute validation statistics for one seed and its optional expansion.

# Arguments

- `target::ResolvedTarget`: resolved target metadata.
- `seed::SeedSelection`: selected ThorAxe seed to validate.
- `expansion`: optional expansion result. Pass `nothing` or `missing` when the
  run stopped after the ThorAxe MSA stage.
- `workdir::AbstractString`: Iduna work directory.

# Keywords

- `overwrite::Bool = false`: reuse package-owned validation outputs when their
  run identity still matches. When `true`, Iduna rebuilds those outputs.

# Returns

- A [`ValidationResult`](@ref) with MSA statistics and optional UniProt comparison
  details.
"""
function validate_results(target::ResolvedTarget,
        seed::SeedSelection,
        expansion::_MaybeExpansion,
        workdir::AbstractString;
        overwrite::Bool = false)
    paths = _validation_input_paths(target, seed, expansion, workdir)
    outputs = _validation_outputs(workdir, seed, paths)
    identity = _validation_identity(target, seed, expansion, paths)
    cache = overwrite ? (; reusable = false, status = :stale, warning = nothing) :
            _classify_stage_state(_validation_dir(workdir, seed), identity, outputs;
        stage_label = "validation")
    if cache.reusable
        @info "Reusing cached validation." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid=seed.pid
        cached = _cached_validation_result(outputs, workdir, seed, expansion)
        _write_validation_state(workdir, seed, :done, identity, outputs;
            action = :reuse,
            warnings = cached.warnings)
        return cached
    end

    cache.warning === nothing ||
        @warn String(cache.warning) validation_dir=_validation_dir(workdir, seed) status=cache.status
    action = cache.status === :missing ? :run : :rebuild
    if overwrite || cache.status !== :missing
        dir = _validation_dir(workdir, seed)
        isdir(dir) && safe_rm(dir, workdir)
    end
    _write_validation_state(workdir, seed, :running, identity, outputs; action)
    try
        stats = _validation_alignment_stats(paths, expansion)
        comparison = _compare_final_query_to_uniprot(
            target, seed, expansion, stats.final_stats, paths.uniprot_path, workdir)
        warnings = comparison.warnings
        stats_path = _write_validation_stats(
            target, seed, stats, comparison, paths, workdir)
        result = _validation_result(stats_path, stats, comparison, warnings)
        _write_validation_state(workdir, seed, :done, identity, outputs; action, warnings)
        return result
    catch err
        _write_validation_state(workdir, seed, :failed, identity, outputs;
            action,
            exception = (;
                type = string(typeof(err)),
                message = sprint(showerror, err)))
        rethrow()
    end
end

end
