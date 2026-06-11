# ThorAxe Result Assembly
# -----------------------

function _selected_seed_msas(seeds::AbstractVector{SeedSelection}, workdir::AbstractString)
    return [read_file(_resolve_artifact_path(seed.stockholm_path, workdir),
                Stockholm; keepinserts = true)
            for seed in seeds]
end

function _selected_translation_warnings(target::ResolvedTarget,
        selected_msas::AbstractVector,
        workdir::AbstractString)
    artifact_workdir = target.workdir === nothing ? workdir : target.workdir
    return collect(Iterators.flatten(
        _validate_transcript_translation(target, selected_msa; workdir = artifact_workdir)
    for selected_msa in selected_msas))
end

function _thoraxe_result_warnings(target::ResolvedTarget,
        seeds::AbstractVector{SeedSelection},
        workdir::AbstractString,
        input_dir::AbstractString,
        summary_path::AbstractString,
        species_filter,
        biomart_filter)
    selected_msas = _selected_seed_msas(seeds, workdir)
    return unique(vcat(species_filter.warnings,
        biomart_filter.warnings,
        _biomart_transcript_query_warnings(input_dir, _thoraxe_logs_dir(workdir)),
        _candidate_summary_warnings(summary_path),
        _selected_translation_warnings(target, selected_msas, workdir)))
end

function _thoraxe_msa_result(input_dir::AbstractString,
        summary_path::AbstractString,
        workdir::AbstractString,
        seeds::AbstractVector{SeedSelection},
        artifacts::AbstractVector,
        warnings::Vector{String};
        pid_sample_count::Integer,
        pid_sample_fraction::Real,
        pid_sample_seed,
        sampling_strategy::Symbol)
    status = isempty(warnings) ? :ok : :warn
    return ThorAxeMSAResult(;
        input_dir,
        thoraxe_dirs = [artifact.thoraxe_dir for artifact in artifacts],
        msa_dir = _thoraxe_msa_dir(workdir),
        baseline_fastas = [String(artifact.seed_fasta) for artifact in artifacts],
        baseline_stockholms = [artifact.seed_path for artifact in artifacts],
        sequence_fastas = [artifact.sequence_fasta for artifact in artifacts],
        species_files = [artifact.species_file for artifact in artifacts],
        pid_summary = summary_path,
        seeds,
        logs_dir = _thoraxe_logs_dir(workdir),
        pid_sample_count = Int(pid_sample_count),
        pid_sample_fraction = Float64(pid_sample_fraction),
        pid_sample_seed,
        sampling_strategy,
        warnings,
        status
    )
end

function _cached_thoraxe_msa_result(input_dir::AbstractString,
        summary_path::AbstractString,
        target::ResolvedTarget,
        workdir::AbstractString,
        cached,
        species_filter,
        biomart_filter;
        pid_sample_count::Integer,
        pid_sample_fraction::Real)
    seeds = cached.seeds
    warnings = _thoraxe_result_warnings(
        target, seeds, workdir, input_dir, summary_path, species_filter, biomart_filter)
    artifacts = [_seed_artifacts(seed, workdir) for seed in seeds]
    return _thoraxe_msa_result(
        input_dir, summary_path, workdir, seeds, artifacts, warnings;
        pid_sample_count,
        pid_sample_fraction,
        pid_sample_seed = cached.sample_seed,
        sampling_strategy = cached.sampling_strategy)
end

# PID Sample Scoring
# ------------------

function _generate_validated_pid_candidate(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid::Real,
        pid_order::Integer,
        specieslist::Union{Nothing, AbstractString};
        sample_count::Integer,
        overwrite::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe)
    candidate_details = (;
        pid,
        pid_order,
        sample_count,
        _specieslist_log_details(specieslist)...)
    @debug "Preparing ThorAxe PID candidate." candidate_details...
    candidate = merge(
        _generate_pid_candidate(target, input_dir, workdir, pid, specieslist;
            overwrite, thoraxe_fn),
        (; workdir))
    validation = _candidate_msa0_validation(target, candidate.msa, pid, workdir)
    if !validation.eligible
        @info "Skipping ineligible ThorAxe PID candidate." pid issue=validation.issue
    end
    return (; pid = Float64(pid), pid_order = Int(pid_order), candidate, validation)
end

function _prepare_candidate_species_samples!(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        records::AbstractVector,
        effective_specieslist::Union{Nothing, AbstractString},
        sampling_strategy::Symbol;
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        overwrite::Bool = false)
    Int(sample_count) == 0 && return nothing
    eligible_records = [record for record in records if record.validation.eligible]
    isempty(eligible_records) && return nothing

    if sampling_strategy === :independent
        @info "Preparing independent ThorAxe PID species samples." sample_fraction
        for record in eligible_records
            # Independent sampling lets each PID candidate use its own available species.
            _ensure_pid_candidate_samples(workdir, record.pid, record.candidate.msa,
                record.candidate.species;
                sample_count,
                sample_fraction,
                sample_seed,
                overwrite,
                gene_id = target.ensembl_gene_id,
                transcript_id = target.transcript_id)
        end
        return nothing
    end

    universe = _shared_sampling_universe(sampling_strategy, target,
        eligible_records, effective_specieslist)
    isempty(universe.species) &&
        error("No species are available for sampling_strategy=$(sampling_strategy).")
    if sampling_strategy === :common
        n_common_species = length(universe.species)
        @info "Preparing common ThorAxe PID species samples." n_common_species sample_fraction
        if n_common_species <= _LOW_COMMON_SPECIES_THRESHOLD
            recommendation = effective_specieslist === nothing ?
                             "consider sampling_strategy=:input or :independent." :
                             "consider sampling_strategy=:input or :independent, or revising the provided specieslist."
            @warn "Common ThorAxe PID species list is small; $(recommendation)" n_common_species warning_threshold=_LOW_COMMON_SPECIES_THRESHOLD sample_fraction
        end
    else
        @info "Preparing input ThorAxe PID species samples." n_species=length(universe.species) sample_fraction
    end
    # Shared strategies write one sampled species list, then point each PID to it.
    _write_shared_species_samples(workdir, universe.species, universe.reference_species;
        sample_count,
        sample_fraction,
        sample_seed,
        overwrite)
    _link_pid_candidate_samples(workdir, [record.pid for record in eligible_records];
        sample_count,
        overwrite)
    return nothing
end

function _pid_sample_spec(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        candidate,
        pid::Real,
        pid_order::Integer,
        sample_idx::Integer,
        n_sequences_reference::Integer;
        thoraxe_fn::Function = ThorAxe.thoraxe)
    sample_paths = _pid_sample_paths(workdir, pid, sample_idx)
    species_file = sample_paths.species_file
    isfile(species_file) || return nothing
    return (;
        gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        pid = Float64(pid),
        pid_order,
        sample = sample_idx,
        sample_label = _candidate_sample_label(sample_idx),
        sample_msa_fasta = sample_paths.fasta_path,
        species_file,
        n_sequences_reference,
        n_sequences_sample = nothing,
        build_sample_msa = (spec;
            overwrite = false) -> begin
            _run_thoraxe_pid_msa(
                target, input_dir, workdir, pid, spec.species_file, sample_idx;
                overwrite, thoraxe_fn)
            return spec.sample_msa_fasta
        end)
end

function _score_pid_samples(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        candidate,
        pid::Real,
        pid_order::Integer;
        sample_count::Integer,
        overwrite::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe,
        score_fn::Function = EPLI.hhsuite_identity_score,
        normalization_fn::Function = EPLI.comparable_positions_normalization,
        progress_output::IO = stderr,
        progress_enabled::Bool = _terminal_progress_enabled(progress_output))
    sample_total = Int(sample_count)
    n_sequences_reference = nsequences(candidate.msa)
    sample_specs = NamedTuple[]
    for sample_idx in 1:sample_total
        spec = _pid_sample_spec(
            target, input_dir, workdir, candidate, pid, pid_order, sample_idx,
            n_sequences_reference; thoraxe_fn)
        spec === nothing || push!(sample_specs, spec)
    end
    isempty(sample_specs) && error("No PID sample scores were computed for PID $(pid).")
    return EPLI._score_alignment_samples(candidate.fasta_path, sample_specs;
        score_fn,
        normalization_fn,
        scores_path = _pid_scores_path(workdir, pid),
        logs_dir = joinpath(_thoraxe_logs_dir(workdir), format_pid_dir(pid)),
        overwrite,
        progress_desc = "Scoring ThorAxe PID $(format_pid(pid)) samples: ",
        progress_output,
        progress_enabled)
end

function _score_validated_pid_candidate(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        record;
        sample_count::Integer,
        sample_fraction::Real,
        sample_seed::UInt64,
        metadata,
        overwrite::Bool = false,
        thoraxe_fn::Function = ThorAxe.thoraxe,
        score_fn::Function = EPLI.hhsuite_identity_score,
        normalization_fn::Function = EPLI.comparable_positions_normalization)
    candidate = record.candidate
    validation = record.validation
    pid = record.pid
    pid_order = record.pid_order
    if !validation.eligible
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

function _score_pid_candidates(target::ResolvedTarget,
        input_dir::AbstractString,
        workdir::AbstractString,
        pid_thresholds::AbstractVector{<:Real},
        effective_specieslist::Union{Nothing, AbstractString},
        metadata;
        pid_sample_count::Integer,
        pid_sample_fraction::Real,
        sample_seed::UInt64,
        overwrite::Bool,
        thoraxe_fn::Function = ThorAxe.thoraxe,
        score_fn::Function = EPLI.hhsuite_identity_score,
        normalization_fn::Function = EPLI.comparable_positions_normalization)
    sampling_strategy = metadata.sampling_strategy
    candidate_records = NamedTuple[]
    for (pid_order, pid) in enumerate(Float64.(pid_thresholds))
        push!(candidate_records,
            _generate_validated_pid_candidate(
                target, input_dir, workdir, pid, pid_order, effective_specieslist;
                sample_count = Int(pid_sample_count),
                overwrite,
                thoraxe_fn))
    end
    _prepare_candidate_species_samples!(target, input_dir, workdir,
        candidate_records, effective_specieslist, sampling_strategy;
        sample_count = Int(pid_sample_count),
        sample_fraction = Float64(pid_sample_fraction),
        sample_seed,
        overwrite)
    return [_score_validated_pid_candidate(
                target, input_dir, workdir, record;
                sample_count = Int(pid_sample_count),
                sample_fraction = Float64(pid_sample_fraction),
                sample_seed,
                metadata,
                overwrite,
                thoraxe_fn,
                score_fn,
                normalization_fn)
            for record in candidate_records]
end

function _select_scored_candidate_seeds(summary_path::AbstractString,
        pid_sample_count::Integer)
    Int(pid_sample_count) == 0 && return _eligible_candidate_seeds(summary_path)
    return [select_best_seed(summary_path)]
end

# Public API
# ----------

"""
    build_thoraxe_msa(target, workdir; pid_thresholds=DEFAULT_PID_THRESHOLDS,
                      specieslist="ases", cached_thoraxe_input_dir=nothing,
                      overwrite=false, orthology="1:1", specieslist_filter=true,
                      biomart_datasets_filter=true, transcript_query_retries=2,
                      pid_sample_count=45, pid_sample_fraction=0.8,
                      pid_sample_seed=nothing,
                      sampling_strategy=:common) -> ThorAxeMSAResult

Run or reuse the ThorAxe input and MSA stages for one resolved target.

# Arguments

- `target::ResolvedTarget`: resolved target metadata.
- `workdir::AbstractString`: Iduna work directory.

# Keywords

- `pid_thresholds = DEFAULT_PID_THRESHOLDS`: ThorAxe percent identity (PID)
  thresholds to test; currently `[10, 20, 30, 60, 80]`.
- `specieslist::AbstractString = "ases"`: species preset, list, file, or name.
  Accepted presets are `"ases"` for the default curated set used on the Ases
  webserver, and `"all"` or `""` for unrestricted ThorAxe species selection.
  Other values are treated as a species-list file path, a comma-separated species
  list, or one species name.
- `cached_thoraxe_input_dir = nothing`: existing transcript-query bundle to copy
  and reuse. When `nothing`, Iduna runs transcript-query to build the input
  bundle.
- `overwrite::Bool = false`: reuse package-owned stage outputs when their run
  identity still matches. When `true`, Iduna rebuilds those outputs.
- `orthology::AbstractString = "1:1"`: Ensembl homology relationship filter.
  Accepted values are `"1:1"` for one-to-one orthologs, `"1:n"` for one-to-one
  and one-to-many orthologs, and `"m:n"` for one-to-one, one-to-many, and
  many-to-many orthologs.
- `specieslist_filter::Bool = true`: filter requested species by Ensembl
  homology. When `false`, Iduna passes the requested species list through without
  this Ensembl homology filter.
- `biomart_datasets_filter::Bool = true`: keep species with BioMart datasets.
  When `false`, Iduna does not remove species that lack a matching BioMart
  dataset.
- `transcript_query_retries::Integer = 2`: number of transcript-query attempts.
- `pid_sample_count::Integer = 45`: number of species samples per PID candidate.
  `0` disables sampling-based seed selection and carries every eligible PID
  candidate forward.
- `pid_sample_fraction::Real = 0.8`: fraction of non-reference species per sample.
  Must be greater than `0` and at most `1`.
- `pid_sample_seed = nothing`: random seed for reproducible PID sampling. When
  `nothing`, Iduna chooses a random seed and records it in the result.
- `sampling_strategy::Symbol = :common`: how species samples are shared across
  PID candidates. Accepted values are `:common`, which samples one shared species
  universe from species common to all eligible PID candidates; `:independent`,
  which samples separately within each PID candidate; and `:input`, which samples
  from the effective input species list after filtering.
"""
function build_thoraxe_msa(target::ResolvedTarget, workdir::AbstractString;
        pid_thresholds::AbstractVector{<:Real} = DEFAULT_PID_THRESHOLDS,
        specieslist::AbstractString = "ases",
        cached_thoraxe_input_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        orthology::AbstractString = "1:1",
        specieslist_filter::Bool = true,
        biomart_datasets_filter::Bool = true,
        transcript_query_retries::Integer = 2,
        pid_sample_count::Integer = 45,
        pid_sample_fraction::Real = 0.8,
        pid_sample_seed::Union{Nothing, Integer} = nothing,
        sampling_strategy::Symbol = :common)
    _orthology_relationships(orthology)
    _validate_pid_sampling_options(pid_sample_count, pid_sample_fraction)
    sampling_strategy = _validate_sampling_strategy(sampling_strategy)
    isempty(pid_thresholds) && error("pid_thresholds cannot be empty.")
    sample_seed = pid_sample_seed === nothing ? UInt64(rand(UInt32)) :
                  _normalize_pid_sample_seed(pid_sample_seed)
    # Species filters run only when transcript_query will create new input.
    resolved_specieslist = _resolve_specieslist_preset(specieslist)
    _log_specieslist_choice(resolved_specieslist)
    requested_species_summary = _specieslist_log_summary(resolved_specieslist.specieslist)
    filters = _resolve_thoraxe_species_filters(
        target, resolved_specieslist.specieslist, orthology,
        cached_thoraxe_input_dir, specieslist_filter, biomart_datasets_filter)
    effective_species_summary = _specieslist_log_summary(filters.effective_specieslist)
    species_details = (;
        n_requested_species = requested_species_summary.n_species,
        n_effective_species = effective_species_summary.n_species,
        specieslist_filter,
        biomart_datasets_filter,
        orthology,
        cached_input = cached_thoraxe_input_dir)
    @info "Resolved ThorAxe transcript_query species." species_details...
    @debug "Resolved ThorAxe transcript_query species previews." requested_specieslist_preview=requested_species_summary.specieslist_preview effective_specieslist_preview=effective_species_summary.specieslist_preview
    input_dir = _ensure_transcript_query(target, workdir;
        specieslist = filters.effective_specieslist, overwrite,
        cached_input_dir = cached_thoraxe_input_dir,
        max_retries = transcript_query_retries,
        orthology)
    summary_path = joinpath(_thoraxe_msa_dir(workdir), "candidate_summary.csv")
    metadata = _candidate_run_metadata(input_dir, target, pid_thresholds;
        sample_count = Int(pid_sample_count),
        sample_fraction = Float64(pid_sample_fraction),
        sample_seed,
        requested_sample_seed = pid_sample_seed,
        sampling_strategy,
        effective_specieslist = filters.effective_specieslist,
        orthology,
        specieslist_filter,
        biomart_datasets_filter)
    stage_identity = _thoraxe_msa_identity(metadata, workdir)
    stage_cache = _thoraxe_msa_stage_cache(workdir, summary_path, metadata,
        stage_identity; overwrite)
    cached_result = _maybe_cached_thoraxe_msa(input_dir, summary_path, target,
        workdir, metadata, stage_identity, stage_cache, filters.species_filter,
        filters.biomart_filter;
        overwrite,
        pid_sample_count,
        pid_sample_fraction)
    cached_result === nothing || return cached_result

    prepared = _prepare_thoraxe_msa_stage!(workdir, summary_path, stage_identity,
        stage_cache; overwrite)
    return _run_thoraxe_msa_stage_with_failure_state!(_run_thoraxe_msa_stage!,
        target, input_dir, workdir, summary_path, pid_thresholds,
        filters.effective_specieslist, metadata, stage_identity, filters, prepared;
        pid_sample_count,
        pid_sample_fraction,
        sample_seed,
        sampling_strategy)
end
