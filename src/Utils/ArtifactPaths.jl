function _is_path_inside_workdir(rel::AbstractString)
    rel == "." && return true
    isabspath(rel) && return false
    parts = splitpath(rel)
    return !isempty(parts) && first(parts) != ".."
end

_relative_artifact_path(path::Nothing, _workdir::AbstractString) = nothing

function _relative_artifact_path(path::AbstractString, workdir::AbstractString)
    str = String(path)
    isabspath(str) || return str
    rel = relpath(abspath(str), abspath(workdir))
    # Only rewrite paths inside workdir; outside inputs must remain absolute.
    return _is_path_inside_workdir(rel) ? rel : str
end

_resolve_artifact_path(path::Nothing, _workdir::AbstractString) = nothing

function _resolve_artifact_path(path::AbstractString, workdir::AbstractString)
    str = String(path)
    return isabspath(str) ? str : joinpath(workdir, str)
end

function _resolve_artifact_path(result::IdunaResult, path)
    _resolve_artifact_path(path, result.workdir)
end

function _relative_seed_paths(seed::SeedSelection, workdir::AbstractString)
    return SeedSelection(;
        pid = seed.pid,
        epli = seed.epli,
        stockholm_path = _relative_artifact_path(seed.stockholm_path, workdir),
        fasta_path = _relative_artifact_path(seed.fasta_path, workdir),
        s_exon_blocks_tsv = _relative_artifact_path(seed.s_exon_blocks_tsv, workdir),
        summary_path = _relative_artifact_path(seed.summary_path, workdir),
        used_fallback_dir = seed.used_fallback_dir,
        workdir
    )
end

function _relative_target_paths(target::ResolvedTarget, workdir::AbstractString)
    return ResolvedTarget(;
        input_id = target.input_id,
        input_kind = target.input_kind,
        uniprot_id = target.uniprot_id,
        ensembl_gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        ensembl_protein_id = target.ensembl_protein_id,
        species = target.species,
        uniprot_sequence_path = _relative_artifact_path(
            target.uniprot_sequence_path, workdir),
        ensembl_protein_sequence_path = _relative_artifact_path(
            target.ensembl_protein_sequence_path, workdir),
        sequence_validated = target.sequence_validated,
        mapping_confirmed = target.mapping_confirmed,
        workdir,
        warnings = target.warnings
    )
end

function _relative_thoraxe_msa_paths(thoraxe::ThorAxeMSAResult, workdir::AbstractString)
    return ThorAxeMSAResult(;
        input_dir = _relative_artifact_path(thoraxe.input_dir, workdir),
        thoraxe_dirs = [_relative_artifact_path(path, workdir)
                        for path in thoraxe.thoraxe_dirs],
        msa_dir = _relative_artifact_path(thoraxe.msa_dir, workdir),
        baseline_fastas = [_relative_artifact_path(path, workdir)
                           for path in thoraxe.baseline_fastas],
        baseline_stockholms = [_relative_artifact_path(path, workdir)
                               for path in thoraxe.baseline_stockholms],
        sequence_fastas = [_relative_artifact_path(path, workdir)
                           for path in thoraxe.sequence_fastas],
        species_files = [_relative_artifact_path(path, workdir)
                         for path in thoraxe.species_files],
        pid_summary = _relative_artifact_path(thoraxe.pid_summary, workdir),
        seeds = [_relative_seed_paths(seed, workdir) for seed in thoraxe.seeds],
        logs_dir = _relative_artifact_path(thoraxe.logs_dir, workdir),
        pid_sample_count = thoraxe.pid_sample_count,
        pid_sample_fraction = thoraxe.pid_sample_fraction,
        pid_sample_seed = thoraxe.pid_sample_seed,
        sampling_strategy = thoraxe.sampling_strategy,
        warnings = thoraxe.warnings,
        status = thoraxe.status
    )
end

function _relative_expansion_paths(expansion::ExpansionResult, workdir::AbstractString)
    return ExpansionResult(;
        run_dir = _relative_artifact_path(expansion.run_dir, workdir),
        seed_stockholm = _relative_artifact_path(expansion.seed_stockholm, workdir),
        seed_fasta = _relative_artifact_path(expansion.seed_fasta, workdir),
        hits_fasta = _relative_artifact_path(expansion.hits_fasta, workdir),
        full_stockholm = _relative_artifact_path(expansion.full_stockholm, workdir),
        match_stockholm = _relative_artifact_path(expansion.match_stockholm, workdir),
        a3m_path = _relative_artifact_path(expansion.a3m_path, workdir),
        s_exon_blocks_tsv = _relative_artifact_path(expansion.s_exon_blocks_tsv, workdir),
        db_dir = _relative_artifact_path(expansion.db_dir, workdir),
        hmm_dir = _relative_artifact_path(expansion.hmm_dir, workdir),
        logs_dir = _relative_artifact_path(expansion.logs_dir, workdir),
        n_hits = expansion.n_hits,
        n_new_hits = expansion.n_new_hits,
        status = expansion.status,
        workdir
    )
end

_relative_expansion_paths(::Missing, _workdir::AbstractString) = missing

function _relative_validation_paths(validation::ValidationResult, workdir::AbstractString)
    return ValidationResult(;
        stats_path = _relative_artifact_path(validation.stats_path, workdir),
        query_name = validation.query_name,
        query_vs_uniprot_path = _relative_artifact_path(
            validation.query_vs_uniprot_path, workdir),
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
        status = validation.status
    )
end

function _relative_result_paths(result::IdunaResult)
    return IdunaResult(;
        input_id = result.input_id,
        workdir = result.workdir,
        target = _relative_target_paths(result.target, result.workdir),
        thoraxe_msa = _relative_thoraxe_msa_paths(result.thoraxe_msa, result.workdir),
        expansions = [_relative_expansion_paths(expansion, result.workdir)
                      for expansion in result.expansions],
        validations = [_relative_validation_paths(validation, result.workdir)
                       for validation in result.validations],
        stages = result.stages,
        warnings = result.warnings,
        status = result.status
    )
end
