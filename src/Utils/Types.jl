"""
    ResolvedTarget

Resolved identifiers and sequence-validation metadata for one input ID.
"""
Base.@kwdef struct ResolvedTarget
    input_id::String
    input_kind::Symbol
    uniprot_id::Union{Nothing, String} = nothing
    ensembl_gene_id::String
    transcript_id::String
    ensembl_protein_id::Union{Nothing, String} = nothing
    species::Union{Nothing, String} = nothing
    uniprot_sequence_path::Union{Nothing, String} = nothing
    ensembl_protein_sequence_path::Union{Nothing, String} = nothing
    sequence_validated::Union{Nothing, Bool} = nothing
    mapping_confirmed::Union{Nothing, Bool} = nothing
    warnings::Vector{String} = String[]
end

"""
    SeedSelection

The ThorAxe PID seed chosen for expansion, plus the summary values used to pick it.
"""
Base.@kwdef struct SeedSelection
    pid::Float64
    median_identity::Float64
    mean_identity::Float64
    stockholm_path::String
    fasta_path::Union{Nothing, String} = nothing
    summary_path::String
    used_fallback_dir::Bool = false
end

"""
    ThorAxeMSAResult

Paths and metadata produced by the ThorAxe MSA-building stage.
"""
Base.@kwdef struct ThorAxeMSAResult
    input_dir::String
    thoraxe_dir::String
    msa_dir::String
    baseline_fasta::String
    baseline_stockholm::String
    sequence_fasta::String
    species_file::String
    pid_summary::String
    best_seed::SeedSelection
    logs_dir::String
    warnings::Vector{String} = String[]
    status::Symbol = :ok
end

"""
    ExpansionResult

Paths and hit counts produced by the MMseqs2/HMMER expansion stage.
"""
Base.@kwdef struct ExpansionResult
    run_dir::String
    seed_stockholm::String
    seed_fasta::Union{Nothing, String} = nothing
    hits_fasta::String
    full_stockholm::String
    match_stockholm::String
    a3m_path::String
    db_dir::String
    hmm_dir::String
    logs_dir::String
    n_hits::Int = 0
    n_new_hits::Int = 0
    status::Symbol = :ok
end

"""
    ValidationResult

Seed and expanded-MSA statistics plus optional query-vs-UniProt checks.
"""
Base.@kwdef struct ValidationResult
    stats_path::String
    query_name::Union{Nothing, String} = nothing
    query_vs_uniprot_path::Union{Nothing, String} = nothing
    seed_nseq::Union{Nothing, Int} = nothing
    seed_ncol::Union{Nothing, Int} = nothing
    seed_clusters62::Union{Nothing, Int} = nothing
    seed_neff80::Union{Nothing, Float64} = nothing
    expanded_nseq::Union{Nothing, Int} = nothing
    expanded_ncol::Union{Nothing, Int} = nothing
    expanded_clusters62::Union{Nothing, Int} = nothing
    expanded_neff80::Union{Nothing, Float64} = nothing
    aln_identical::Union{Nothing, Bool} = nothing
    aln_mismatches::Union{Nothing, Int} = nothing
    aln_insertions::Union{Nothing, Int} = nothing
    aln_deletions::Union{Nothing, Int} = nothing
    warnings::Vector{String} = String[]
    status::Symbol = :ok
end

"""
    IdunaResult

Top-level result returned by `Iduna.iduna`. It keeps the full pipeline status
and links to the per-stage result objects. `expansion` is `nothing` when the
run was requested with `no_expansion=true`.
"""
Base.@kwdef struct IdunaResult
    input_id::String
    workdir::String
    target::ResolvedTarget
    thoraxe_msa::ThorAxeMSAResult
    expansion::Union{Nothing, ExpansionResult}
    validation::ValidationResult
    warnings::Vector{String} = String[]
    status::Symbol = :ok
end
