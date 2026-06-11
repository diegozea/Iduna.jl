# Stage Summaries
# ---------------

function _stage_summary_from_state_path(state_path::AbstractString, workdir::AbstractString)
    state = try
        JSON.parse(read(state_path, String))
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    state isa AbstractDict || return nothing
    return (;
        stage = get(state, "stage", nothing),
        stage_key = get(state, "stage_key", nothing),
        status = get(state, "status", nothing),
        action = get(state, "action", nothing),
        identity_hash = get(state, "identity_hash", nothing),
        state_path = _relative_artifact_path(state_path, workdir),
        outputs = get(state, "outputs", Dict{String, Any}()),
        warnings = get(state, "warnings", String[]),
        exception = get(state, "exception", nothing)
    )
end

function _push_stage_state_paths!(paths::Vector{String}, root::AbstractString)
    isdir(root) || return paths
    for (dir, _, files) in walkdir(root)
        _STAGE_STATE_FILE in files || continue
        push!(paths, joinpath(dir, _STAGE_STATE_FILE))
    end
    return paths
end

"""
    collect_stage_summaries(workdir; stage_keys=nothing)

Read stage summaries from an Iduna work directory.

# Arguments

- `workdir::AbstractString`: Iduna work directory to scan.

# Keywords

- `stage_keys = nothing`: optional stage keys to keep. When `nothing`, all found
  stage summaries are returned.
"""
function collect_stage_summaries(workdir::AbstractString; stage_keys = nothing)
    paths = String[]
    _push_stage_state_paths!(paths, joinpath(workdir, ".iduna", "stages"))
    _push_stage_state_paths!(paths, joinpath(workdir, "expansion"))
    _push_stage_state_paths!(paths, joinpath(workdir, "validation"))
    stage_key_set = stage_keys === nothing ? nothing : Set(String.(stage_keys))
    summaries = Any[]
    for path in sort!(unique(paths))
        summary = _stage_summary_from_state_path(path, workdir)
        summary === nothing && continue
        if stage_key_set !== nothing
            summary.stage_key === nothing && continue
            String(summary.stage_key) in stage_key_set || continue
        end
        push!(summaries, summary)
    end
    return summaries
end

# Result Summaries
# ----------------

function _seed_summary(seed::SeedSelection)
    return (;
        pid = seed.pid,
        epli = seed.epli === missing ? nothing : seed.epli,
        stockholm_path = seed.stockholm_path,
        fasta_path = seed.fasta_path,
        s_exon_blocks_tsv = seed.s_exon_blocks_tsv,
        summary_path = seed.summary_path,
        used_fallback_dir = seed.used_fallback_dir
    )
end

function _expansion_summary(expansion::ExpansionResult)
    return (;
        run_dir = expansion.run_dir,
        seed_stockholm = expansion.seed_stockholm,
        seed_fasta = expansion.seed_fasta,
        match_stockholm = expansion.match_stockholm,
        full_stockholm = expansion.full_stockholm,
        a3m_path = expansion.a3m_path,
        s_exon_blocks_tsv = expansion.s_exon_blocks_tsv,
        hits_fasta = expansion.hits_fasta,
        db_dir = expansion.db_dir,
        hmm_dir = expansion.hmm_dir,
        logs_dir = expansion.logs_dir,
        n_hits = expansion.n_hits,
        n_new_hits = expansion.n_new_hits,
        status = String(expansion.status)
    )
end

_expansion_summary(::Missing) = nothing

function _validation_summary(validation::ValidationResult)
    return (;
        stats_path = validation.stats_path,
        query_name = validation.query_name,
        query_vs_uniprot_path = validation.query_vs_uniprot_path,
        seed_nseq = validation.seed_nseq,
        seed_ncol = validation.seed_ncol,
        seed_clusters62 = validation.seed_clusters62,
        seed_neff80 = validation.seed_neff80,
        expanded_nseq = validation.expanded_nseq,
        expanded_ncol = validation.expanded_ncol,
        expanded_clusters62 = validation.expanded_clusters62,
        expanded_neff80 = validation.expanded_neff80,
        aln_identical = validation.aln_identical,
        aln_mismatches = validation.aln_mismatches,
        aln_insertions = validation.aln_insertions,
        aln_deletions = validation.aln_deletions,
        warnings = validation.warnings,
        status = String(validation.status)
    )
end

"""
    result_summary(result) -> NamedTuple

Return a compact, JSON-friendly summary of an [`IdunaResult`](@ref).

# Arguments

- `result::IdunaResult`: result to summarize.
"""
function result_summary(result::IdunaResult)
    result = _relative_result_paths(result)
    return (;
        input_id = result.input_id,
        status = String(result.status),
        warnings = result.warnings,
        stages = result.stages,
        target = (;
            input_kind = String(result.target.input_kind),
            uniprot_id = result.target.uniprot_id,
            ensembl_gene_id = result.target.ensembl_gene_id,
            transcript_id = result.target.transcript_id,
            ensembl_protein_id = result.target.ensembl_protein_id,
            species = result.target.species,
            uniprot_sequence_path = result.target.uniprot_sequence_path,
            ensembl_protein_sequence_path = result.target.ensembl_protein_sequence_path,
            sequence_validated = result.target.sequence_validated,
            mapping_confirmed = result.target.mapping_confirmed
        ),
        thoraxe_msa = (;
            baseline_fastas = result.thoraxe_msa.baseline_fastas,
            baseline_stockholms = result.thoraxe_msa.baseline_stockholms,
            pid_summary = result.thoraxe_msa.pid_summary,
            seeds = _seed_summary.(result.thoraxe_msa.seeds),
            selected_pids = [seed.pid for seed in result.thoraxe_msa.seeds],
            seed_stockholms = [seed.stockholm_path for seed in result.thoraxe_msa.seeds],
            s_exon_blocks_tsvs = [seed.s_exon_blocks_tsv
                                  for seed in result.thoraxe_msa.seeds],
            pid_sample_count = result.thoraxe_msa.pid_sample_count,
            pid_sample_fraction = result.thoraxe_msa.pid_sample_fraction,
            pid_sample_seed = result.thoraxe_msa.pid_sample_seed,
            sampling_strategy = String(result.thoraxe_msa.sampling_strategy),
            warnings = result.thoraxe_msa.warnings,
            status = String(result.thoraxe_msa.status)
        ),
        expansions = _expansion_summary.(result.expansions),
        validations = _validation_summary.(result.validations)
    )
end
