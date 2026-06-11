# Candidate Scoring
# -----------------

function _candidate_msa0_validation(target::ResolvedTarget,
        msa::AbstractMultipleSequenceAlignment,
        pid::Real,
        workdir::AbstractString)
    try
        validation_warnings = _validate_transcript_translation(
            target, msa; workdir = target.workdir === nothing ? workdir : target.workdir)
        issue = isempty(validation_warnings) ? missing : join(validation_warnings, " | ")
        warnings = ["PID $(format_pid(pid)) candidate retained with warning: $(warning)"
                    for warning in validation_warnings]
        status = isempty(validation_warnings) ? "ok" : "warning"
        return (; eligible = true, status, issue, warnings)
    catch err
        err isa InterruptException && rethrow()
        issue = sprint(showerror, err)
        warning = "PID $(format_pid(pid)) excluded from seed selection: $(issue)"
        return (; eligible = false, status = "invalid_msa0", issue, warnings = [warning])
    end
end

_candidate_summary_optional(value) = value === nothing ? missing : value

function _candidate_summary_path(path, workdir::AbstractString)
    (path === nothing || ismissing(path)) && return missing
    return _relative_artifact_path(String(path), workdir)
end

function _candidate_summary_row(target::ResolvedTarget,
        candidate,
        pid::Real,
        pid_order::Integer,
        validation;
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        metadata,
        epli = missing,
        n_samples::Integer = 0)
    return (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        pid = Float64(pid),
        pid_order,
        eligible = Bool(validation.eligible),
        selected = false,
        msa0_status = String(validation.status),
        msa0_issue = validation.issue,
        epli,
        n_samples = Int(n_samples),
        n_sequences_msa0 = nsequences(candidate.msa),
        pid_sample_count = Int(sample_count),
        pid_sample_fraction = Float64(sample_fraction),
        pid_sample_seed = sample_seed,
        sampling_strategy = String(metadata.sampling_strategy),
        pid_thresholds_key = metadata.pid_thresholds_key,
        effective_specieslist = _candidate_summary_optional(metadata.effective_specieslist),
        orthology = metadata.orthology,
        specieslist_filter = metadata.specieslist_filter,
        biomart_datasets_filter = metadata.biomart_datasets_filter,
        transcript_query_fingerprint = _candidate_summary_optional(
            metadata.transcript_query_fingerprint),
        selection_mode = metadata.selection_mode,
        fasta_path = _candidate_summary_path(candidate.fasta_path, candidate.workdir),
        stockholm_path = _candidate_summary_path(candidate.stockholm_path, candidate.workdir),
        sequence_fasta = _candidate_summary_path(candidate.sequence_fasta, candidate.workdir),
        species_file = _candidate_summary_path(candidate.species_file, candidate.workdir),
        scores_path = _candidate_summary_path(_pid_scores_path(candidate.workdir, pid),
            candidate.workdir)
    )
end

function _score_pid_candidate(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        pid_order::Integer,
        specieslist::Union{Nothing, AbstractString};
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        metadata,
        overwrite::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe,
        score_fn::Function = EPLI.hhsuite_identity_score,
        normalization_fn::Function = EPLI.comparable_positions_normalization)
    @debug "Scoring ThorAxe PID candidate." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id pid pid_order sample_count
    candidate = merge(
        _generate_pid_candidate(target, input_dir, workdir, pid, specieslist;
            overwrite, thoraxe_fn),
        (; workdir))
    validation = _candidate_msa0_validation(target, candidate.msa, pid, workdir)
    if !validation.eligible
        @info "Skipping ineligible ThorAxe PID candidate." pid issue=validation.issue
        return _candidate_summary_row(
            target, candidate, pid, pid_order, validation;
            sample_count, sample_fraction, sample_seed, metadata)
    end

    if sample_count == 0
        @info "Keeping ThorAxe PID candidate without sample scoring." pid
        return _candidate_summary_row(
            target, candidate, pid, pid_order, validation;
            sample_count, sample_fraction, sample_seed, metadata)
    end

    _ensure_pid_candidate_samples(workdir, pid, candidate.msa, candidate.species;
        sample_count,
        sample_fraction,
        sample_seed,
        overwrite,
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id)

    scored = _score_pid_samples(
        target, input_dir, workdir, candidate, pid, pid_order;
        sample_count,
        overwrite,
        thoraxe_fn,
        score_fn,
        normalization_fn)
    return _candidate_summary_row(
        target, candidate, pid, pid_order, validation;
        epli = scored.epli,
        n_samples = length(scored.rows),
        sample_count,
        sample_fraction,
        sample_seed,
        metadata)
end

function _summarize_candidate_scores(rows::AbstractVector{<:NamedTuple},
        path::AbstractString)
    df = DataFrame(rows)
    if :pid_order in propertynames(df)
        sort!(df, [:pid_order])
    end
    CSV.write(path, df)
    return path
end

_truthy(value) = value === true || value == 1 || lowercase(string(value)) == "true"
