# Candidate Summary Cache
# -----------------------

const _CANDIDATE_SUMMARY_STRING_TYPES = Dict(
    :pid_thresholds_key => String,
    :effective_specieslist => String,
    :orthology => String,
    :transcript_query_fingerprint => String,
    :selection_mode => String,
    :sampling_strategy => String,
    :msa0_status => String,
    :msa0_issue => String,
    :stockholm_path => String,
    :fasta_path => String,
    :sequence_fasta => String,
    :species_file => String,
    :scores_path => String
)

function _candidate_summary_file(path::AbstractString)
    CSV.File(path; types = _CANDIDATE_SUMMARY_STRING_TYPES, validate = false)
end

function _candidate_summary_dataframe(path::AbstractString)
    DataFrame(_candidate_summary_file(path))
end

function _pid_thresholds_key(pid_thresholds::AbstractVector{<:Real})
    return join(repr.(Float64.(pid_thresholds)), ",")
end

function _candidate_run_metadata(input_dir::AbstractString,
        target::ResolvedTarget,
        pid_thresholds::AbstractVector{<:Real};
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        requested_sample_seed::Union{Nothing, Integer},
        sampling_strategy::Symbol,
        effective_specieslist::Union{Nothing, AbstractString},
        orthology::AbstractString,
        specieslist_filter::Bool,
        biomart_datasets_filter::Bool)
    return (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        pid_thresholds_key = _pid_thresholds_key(pid_thresholds),
        pid_sample_count = Int(sample_count),
        pid_sample_fraction = Float64(sample_fraction),
        pid_sample_seed = sample_seed,
        requested_pid_sample_seed = requested_sample_seed,
        sampling_strategy = sampling_strategy,
        effective_specieslist = _normalized_specieslist(effective_specieslist),
        orthology = String(orthology),
        specieslist_filter = Bool(specieslist_filter),
        biomart_datasets_filter = Bool(biomart_datasets_filter),
        transcript_query_fingerprint = _bundle_fingerprint(input_dir),
        selection_mode = sample_count == 0 ? "all_candidates" : "sampled_selection"
    )
end

function _thoraxe_input_identity_hash(workdir::AbstractString)
    state = _read_stage_state(_thoraxe_input_stage_dir(workdir))
    state isa AbstractDict || return nothing
    value = get(state, "identity_hash", nothing)
    return value === nothing ? nothing : String(value)
end

function _canonical_thoraxe_msa_metadata(metadata)
    if metadata.requested_pid_sample_seed === nothing
        # Auto-generated seeds should not force a rebuild of already scored candidates.
        return merge(metadata, (; pid_sample_seed = nothing))
    end
    return metadata
end

function _thoraxe_msa_identity(metadata, workdir::AbstractString)
    return merge(_canonical_thoraxe_msa_metadata(metadata),
        (; thoraxe_input_identity_hash = _thoraxe_input_identity_hash(workdir)))
end

function _thoraxe_msa_outputs(summary_path::AbstractString,
        seeds::AbstractVector{SeedSelection},
        workdir::AbstractString)
    outputs = Dict{String, Any}("pid_summary" => summary_path)
    outputs["seed_stockholms"] = [_resolve_artifact_path(seed.stockholm_path, workdir)
                                  for seed in seeds]
    outputs["seed_fastas"] = [_resolve_artifact_path(seed.fasta_path, workdir)
                              for seed in seeds if seed.fasta_path !== nothing]
    return outputs
end

function _write_thoraxe_msa_state(workdir::AbstractString,
        summary_path::AbstractString,
        seeds::AbstractVector{SeedSelection},
        status::Symbol,
        identity;
        action,
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing)
    return _write_stage_state(_thoraxe_msa_stage_dir(workdir);
        stage = "thoraxe_msa",
        stage_key = "thoraxe_msa",
        status,
        identity,
        outputs = _thoraxe_msa_outputs(summary_path, seeds, workdir),
        action,
        warnings,
        exception,
        workdir)
end

function _thoraxe_msa_stage_cache(workdir::AbstractString,
        summary_path::AbstractString,
        metadata,
        stage_identity;
        overwrite::Bool)
    cache = overwrite ? (; reusable = false, status = :stale, warning = nothing) :
            _classify_stage_state(_thoraxe_msa_stage_dir(workdir), stage_identity,
        (; pid_summary = summary_path); stage_label = "ThorAxe MSA")
    has_manifest = isfile(joinpath(_thoraxe_msa_stage_dir(workdir), "stage_state.json"))
    summary_matches = !overwrite && _has_matching_candidate_summary(summary_path, metadata)
    return (; cache, has_manifest, summary_matches)
end

function _maybe_cached_thoraxe_msa(input_dir::AbstractString,
        summary_path::AbstractString,
        target::ResolvedTarget,
        workdir::AbstractString,
        metadata,
        stage_identity,
        stage_cache,
        species_filter,
        biomart_filter;
        overwrite::Bool,
        pid_sample_count::Integer,
        pid_sample_fraction::Real)
    can_try_cache = stage_cache.cache.reusable ||
                    (!overwrite && stage_cache.summary_matches)
    can_try_cache || return nothing
    cached = _cached_selected_seeds(summary_path, workdir, metadata)
    cached === nothing && return nothing
    @info "Reusing cached ThorAxe MSA candidates." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id n_seeds=length(cached.seeds) summary_path adopted_legacy=(!stage_cache.cache.reusable&&!stage_cache.has_manifest)
    result = _cached_thoraxe_msa_result(
        input_dir, summary_path, target, workdir, cached,
        species_filter, biomart_filter;
        pid_sample_count,
        pid_sample_fraction)
    _write_thoraxe_msa_state(workdir, summary_path, result.seeds, :done,
        stage_identity; action = :reuse, warnings = result.warnings)
    return result
end

function _prepare_thoraxe_msa_stage!(workdir::AbstractString,
        summary_path::AbstractString,
        stage_identity,
        stage_cache;
        overwrite::Bool)
    if stage_cache.cache.reusable
        @warn "ThorAxe MSA manifest matched, but selected seed artifacts were incomplete; rebuilding." summary_path
    elseif stage_cache.cache.warning !== nothing
        @warn String(stage_cache.cache.warning) summary_path status=stage_cache.cache.status
    end
    action = _stage_action(stage_cache.cache)
    if overwrite || (stage_cache.cache.status !== :missing && stage_cache.has_manifest)
        isdir(_thoraxe_msa_dir(workdir)) && safe_rm(_thoraxe_msa_dir(workdir), workdir)
    end
    # A legacy summary can still be reused when all required files and metadata match.
    local_artifacts_are_current = !overwrite && !stage_cache.has_manifest &&
                                  stage_cache.summary_matches
    force_pid_rerun = overwrite || !local_artifacts_are_current
    _write_thoraxe_msa_state(workdir, summary_path, SeedSelection[], :running,
        stage_identity; action)
    return (; action, local_artifacts_are_current, force_pid_rerun)
end

function _run_thoraxe_msa_stage!(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        summary_path::AbstractString,
        pid_thresholds::AbstractVector{<:Real},
        effective_specieslist::Union{Nothing, AbstractString},
        metadata,
        stage_identity,
        filters,
        prepared;
        pid_sample_count::Integer,
        pid_sample_fraction::Real,
        sample_seed::UInt64,
        sampling_strategy::Symbol)
    build_details = (;
        n_pids = length(pid_thresholds),
        pid_thresholds_preview = _pid_thresholds_preview(pid_thresholds),
        pid_sample_count,
        sampling_strategy)
    @info "Building ThorAxe MSA candidates." build_details...
    score_rows = _score_pid_candidates(target, input_dir, workdir, pid_thresholds,
        effective_specieslist, metadata;
        pid_sample_count,
        pid_sample_fraction,
        sample_seed,
        overwrite = prepared.force_pid_rerun)
    _summarize_candidate_scores(score_rows, summary_path)
    @info "Selecting ThorAxe seed candidates." summary_path pid_sample_count
    seeds = _select_scored_candidate_seeds(summary_path, pid_sample_count)
    _mark_selected_candidates!(summary_path, seeds)
    warnings = _thoraxe_result_warnings(target, seeds, workdir, input_dir, summary_path,
        filters.species_filter, filters.biomart_filter)
    artifacts = [_seed_artifacts(seed, workdir) for seed in seeds]
    _write_thoraxe_msa_state(workdir, summary_path, seeds, :done,
        stage_identity; action = prepared.action, warnings)
    return _thoraxe_msa_result(
        input_dir, summary_path, workdir, seeds, artifacts, warnings;
        pid_sample_count,
        pid_sample_fraction,
        pid_sample_seed = sample_seed,
        sampling_strategy)
end

function _run_thoraxe_msa_stage_with_failure_state!(stage_runner::Function,
        target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        summary_path::AbstractString,
        pid_thresholds::AbstractVector{<:Real},
        effective_specieslist::Union{Nothing, AbstractString},
        metadata,
        stage_identity,
        filters,
        prepared;
        pid_sample_count::Integer,
        pid_sample_fraction::Real,
        sample_seed::UInt64,
        sampling_strategy::Symbol)
    try
        return stage_runner(target, input_dir, workdir, summary_path, pid_thresholds,
            effective_specieslist, metadata, stage_identity, filters, prepared;
            pid_sample_count,
            pid_sample_fraction,
            sample_seed,
            sampling_strategy)
    catch err
        _write_thoraxe_msa_state(workdir, summary_path, SeedSelection[], :failed,
            stage_identity;
            action = prepared.action,
            exception = _exception_summary(err))
        rethrow()
    end
end

function _candidate_summary_matches(df::DataFrame, metadata)
    isempty(df) && return false
    required = (
        :gene_id,
        :transcript_id,
        :pid,
        :pid_thresholds_key,
        :pid_sample_count,
        :pid_sample_fraction,
        :pid_sample_seed,
        :sampling_strategy,
        :effective_specieslist,
        :orthology,
        :specieslist_filter,
        :biomart_datasets_filter,
        :transcript_query_fingerprint,
        :selection_mode
    )
    all(name -> name in propertynames(df), required) || return false
    actual_pids = Float64.(df.pid)
    expected_pids = parse.(Float64, split(metadata.pid_thresholds_key, ','))
    if length(actual_pids) != length(expected_pids) ||
       Set(actual_pids) != Set(expected_pids)
        return false
    end
    return all(eachrow(df)) do row
        checks = (
            string(row.gene_id) == metadata.gene_id,
            string(row.transcript_id) == metadata.transcript_id,
            string(row.pid_thresholds_key) == metadata.pid_thresholds_key,
            Int(row.pid_sample_count) == metadata.pid_sample_count,
            Float64(row.pid_sample_fraction) == metadata.pid_sample_fraction,
            metadata.requested_pid_sample_seed === nothing ||
            UInt64(row.pid_sample_seed) == metadata.pid_sample_seed,
            string(row.sampling_strategy) == String(metadata.sampling_strategy),
            _metadata_value_matches(row.effective_specieslist,
                metadata.effective_specieslist),
            string(row.orthology) == metadata.orthology,
            _truthy(row.specieslist_filter) == metadata.specieslist_filter,
            _truthy(row.biomart_datasets_filter) == metadata.biomart_datasets_filter,
            string(row.transcript_query_fingerprint) ==
            string(metadata.transcript_query_fingerprint),
            string(row.selection_mode) == metadata.selection_mode
        )
        all(checks)
    end
end

function _row_seed(row, summary_path::AbstractString)
    epli = :epli in propertynames(row) && !ismissing(row.epli) ?
           Float64(row.epli) : missing
    workdir = abspath(_candidate_summary_workdir(summary_path))
    return SeedSelection(;
        pid = Float64(row.pid),
        epli,
        stockholm_path = String(row.stockholm_path),
        fasta_path = row.fasta_path === missing ? nothing : String(row.fasta_path),
        s_exon_blocks_tsv = s_exon_blocks_path(String(row.stockholm_path)),
        summary_path,
        workdir
    )
end
