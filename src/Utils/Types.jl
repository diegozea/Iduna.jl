using AutoPrettyPrinting: @def_pprint

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
    workdir::Union{Nothing, String} = nothing
    warnings::Vector{String} = String[]
end

"""
    SeedSelection

The ThorAxe PID seed chosen for expansion, plus the summary values used to pick it.
"""
Base.@kwdef struct SeedSelection
    pid::Float64
    median_identity::Union{Missing, Float64}
    mean_identity::Union{Missing, Float64}
    stockholm_path::String
    fasta_path::Union{Nothing, String} = nothing
    s_exon_blocks_tsv::Union{Nothing, String} = nothing
    summary_path::String
    used_fallback_dir::Bool = false
    workdir::Union{Nothing, String} = nothing
end

"""
    ThorAxeMSAResult

Paths and metadata produced by the ThorAxe MSA-building stage.
"""
Base.@kwdef struct ThorAxeMSAResult
    input_dir::String
    thoraxe_dirs::Vector{String}
    msa_dir::String
    baseline_fastas::Vector{String}
    baseline_stockholms::Vector{String}
    sequence_fastas::Vector{String}
    species_files::Vector{String}
    pid_summary::String
    seeds::Vector{SeedSelection}
    logs_dir::String
    pid_sample_count::Int = 0
    pid_sample_fraction::Float64 = 1.0
    pid_sample_seed::UInt64 = UInt64(0)
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
    s_exon_blocks_tsv::Union{Nothing, String} = nothing
    db_dir::String
    hmm_dir::String
    logs_dir::String
    n_hits::Int = 0
    n_new_hits::Int = 0
    status::Symbol = :ok
    workdir::Union{Nothing, String} = nothing
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
and links to the per-stage result objects. `expansions` is empty when the run
was requested with `no_expansion=true`.
"""
Base.@kwdef struct IdunaResult
    input_id::String
    workdir::String
    target::ResolvedTarget
    thoraxe_msa::ThorAxeMSAResult
    expansions::Vector{ExpansionResult}
    validations::Vector{ValidationResult}
    stages::Vector{Any} = Any[]
    warnings::Vector{String} = String[]
    status::Symbol = :ok
end

@def_pprint mime_types="text/plain" base_show=true ResolvedTarget
@def_pprint mime_types="text/plain" base_show=true SeedSelection
@def_pprint mime_types="text/plain" base_show=true ThorAxeMSAResult
@def_pprint mime_types="text/plain" base_show=true ExpansionResult
@def_pprint mime_types="text/plain" base_show=true ValidationResult
@def_pprint mime_types="text/plain" base_show=true IdunaResult
