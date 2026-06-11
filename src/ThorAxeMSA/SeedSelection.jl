# Seed Selection
# --------------

"""
    select_best_seed(summary_path) -> SeedSelection

Choose the best eligible percent identity (PID) candidate from a ThorAxe
candidate summary table.

Candidates are ranked by EPLI, candidate MSA size, and then the original PID
order.

# Arguments

- `summary_path::AbstractString`: ThorAxe candidate summary CSV path.
"""
function select_best_seed(summary_path::AbstractString)
    df = _candidate_summary_dataframe(summary_path)
    isempty(df) && error("Cannot select a seed from an empty candidate summary.")
    if :eligible in propertynames(df)
        df = df[[!ismissing(value) && _truthy(value) for value in df.eligible], :]
    end
    :epli in propertynames(df) || error("Candidate summary has no epli column.")
    df = df[[!ismissing(value) for value in df.epli], :]
    isempty(df) && error("No eligible PID candidates are available for seed selection.")
    df.__row_order = 1:nrow(df)
    order_col = :pid_order in propertynames(df) ? :pid_order : :__row_order
    nseq_col = if :n_sequences_msa0 in propertynames(df)
        :n_sequences_msa0
    else
        df.__n_sequences_msa0 = zeros(Int, nrow(df))
        :__n_sequences_msa0
    end
    sort!(df, [:epli, nseq_col, order_col]; rev = [true, true, false])
    row = first(eachrow(df))
    return _row_seed(row, summary_path)
end

function _candidate_summary_workdir(summary_path::AbstractString)
    summary_dir = dirname(summary_path)
    return basename(summary_dir) == "thoraxe_msa" ? dirname(summary_dir) : summary_dir
end

function _same_candidate_path(row_path, seed_path, workdir::AbstractString)
    row_path === missing && return false
    row_string = String(row_path)
    seed_string = String(seed_path)
    row_string == seed_string && return true
    return _resolve_artifact_path(row_string, workdir) ==
           _resolve_artifact_path(seed_string, workdir)
end

function _mark_selected_candidates!(summary_path::AbstractString,
        seeds::AbstractVector{SeedSelection})
    df = _candidate_summary_dataframe(summary_path)
    df.selected = falses(nrow(df))
    workdir = _candidate_summary_workdir(summary_path)
    for seed in seeds
        selected_idx = findfirst(eachrow(df)) do row
            Float64(row.pid) == seed.pid &&
                _same_candidate_path(row.stockholm_path, seed.stockholm_path, workdir)
        end
        selected_idx === nothing ||
            (df.selected[selected_idx] = true)
    end
    CSV.write(summary_path, df)
    return summary_path
end

function _mark_selected_candidate!(summary_path::AbstractString, seed::SeedSelection)
    _mark_selected_candidates!(summary_path, [seed])
end

function _selected_candidate_seeds(summary_path::AbstractString)
    df = _candidate_summary_dataframe(summary_path)
    :selected in propertynames(df) || return [select_best_seed(summary_path)]
    selected = df[[!ismissing(value) && _truthy(value) for value in df.selected], :]
    isempty(selected) && return SeedSelection[]
    return [_row_seed(row, summary_path) for row in eachrow(selected)]
end

function _eligible_candidate_seeds(summary_path::AbstractString)
    df = _candidate_summary_dataframe(summary_path)
    :eligible in propertynames(df) || error("Candidate summary has no eligible column.")
    eligible = df[[!ismissing(value) && _truthy(value) for value in df.eligible], :]
    isempty(eligible) && error("No eligible PID candidates are available.")
    sort!(eligible, [:pid_order])
    return [_row_seed(row, summary_path) for row in eachrow(eligible)]
end

function _normalize_pid_sample_seed(seed::Integer)::UInt64
    seed < 0 && error("pid_sample_seed must be non-negative.")
    return UInt64(seed)
end

function _validate_pid_sampling_options(sample_count::Integer, sample_fraction::Real)
    sample_count >= 0 || error("pid_sample_count must be non-negative.")
    0.0 < Float64(sample_fraction) <= 1.0 ||
        error("pid_sample_fraction must be greater than 0 and at most 1.")
    return nothing
end

function _validate_sampling_strategy(sampling_strategy::Symbol)
    sampling_strategy in _SAMPLING_STRATEGIES ||
        error("sampling_strategy must be one of :independent, :common, or :input.")
    return sampling_strategy
end

function _summary_seed_value(df::DataFrame, fallback::UInt64)
    :pid_sample_seed in propertynames(df) || return fallback
    isempty(df) && return fallback
    value = df.pid_sample_seed[1]
    value === missing && return fallback
    return UInt64(value)
end

function _summary_sampling_strategy(df::DataFrame, fallback::Symbol)
    :sampling_strategy in propertynames(df) || return fallback
    isempty(df) && return fallback
    value = df.sampling_strategy[1]
    value === missing && return fallback
    return Symbol(String(value))
end

function _has_current_candidate_summary(df::DataFrame)
    required = Set([
        :gene_id,
        :transcript_id,
        :pid,
        :eligible,
        :selected,
        :msa0_status,
        :msa0_issue,
        :epli,
        :n_sequences_msa0,
        :pid_thresholds_key,
        :sampling_strategy,
        :effective_specieslist,
        :orthology,
        :specieslist_filter,
        :biomart_datasets_filter,
        :transcript_query_fingerprint,
        :selection_mode,
        :sequence_fasta,
        :species_file,
        :scores_path
    ])
    names = Set(propertynames(df))
    return issubset(required, names)
end

function _candidate_summary_warnings(summary_path::AbstractString)
    isfile(summary_path) || return String[]
    df = _candidate_summary_dataframe(summary_path)
    _has_current_candidate_summary(df) || return String[]
    warnings = String[]
    for row in eachrow(df)
        issue = row.msa0_issue
        (ismissing(issue) || isempty(String(issue))) && continue
        pid_label = format_pid(Float64(row.pid))
        if !(!ismissing(row.eligible) && _truthy(row.eligible))
            push!(warnings, "PID $(pid_label) excluded from seed selection: $(issue)")
        elseif String(row.msa0_status) != "ok"
            push!(warnings, "PID $(pid_label) candidate retained with warning: $(issue)")
        end
    end
    return warnings
end

function _seed_artifacts(seed::SeedSelection, workdir::AbstractString)
    paths = _pid_sample_paths(workdir, seed.pid, 0)
    return (;
        seed_path = _resolve_artifact_path(seed.stockholm_path, workdir),
        seed_fasta = seed.fasta_path === nothing ? nothing :
                     _resolve_artifact_path(seed.fasta_path, workdir),
        s_exon_blocks_tsv = seed.s_exon_blocks_tsv === nothing ?
                            s_exon_blocks_path(
            _resolve_artifact_path(seed.stockholm_path, workdir)) :
                            _resolve_artifact_path(seed.s_exon_blocks_tsv, workdir),
        sequence_fasta = paths.sequence_fasta,
        species_file = paths.species_file,
        thoraxe_dir = _pid_sample_thoraxe_dir(workdir, seed.pid, 0)
    )
end

function _selected_artifacts_exist(seed::SeedSelection, workdir::AbstractString)
    artifacts = _seed_artifacts(seed, workdir)
    path_table = joinpath(artifacts.thoraxe_dir, "path_table.csv")
    return isfile(artifacts.seed_path) &&
           artifacts.seed_fasta !== nothing &&
           isfile(artifacts.seed_fasta) &&
           isfile(artifacts.sequence_fasta) &&
           isfile(artifacts.species_file) &&
           isfile(path_table) &&
           _seed_has_s_exon_annotations(artifacts.seed_path)
end

function _cached_selected_seeds(summary_path::AbstractString, workdir::AbstractString,
        metadata)
    isfile(summary_path) || return nothing
    df = _candidate_summary_dataframe(summary_path)
    isempty(df) && return nothing
    _has_current_candidate_summary(df) || return nothing
    _candidate_summary_matches(df, metadata) || return nothing
    seeds = _selected_candidate_seeds(summary_path)
    isempty(seeds) && return nothing
    all(seed -> _selected_artifacts_exist(seed, workdir), seeds) || return nothing
    for seed in seeds
        _ensure_seed_blocks_tsv(_pid_sample_paths(workdir, seed.pid, 0), seed.pid)
    end
    return (;
        seeds,
        sample_seed = _summary_seed_value(df, metadata.pid_sample_seed),
        sampling_strategy = _summary_sampling_strategy(df, metadata.sampling_strategy))
end

function _has_matching_candidate_summary(summary_path::AbstractString, metadata)
    isfile(summary_path) || return false
    df = _candidate_summary_dataframe(summary_path)
    isempty(df) && return false
    _has_current_candidate_summary(df) || return false
    return _candidate_summary_matches(df, metadata)
end

function _resolve_thoraxe_species_filters(target::ResolvedTarget,
        specieslist::Union{Nothing, AbstractString},
        orthology::AbstractString,
        cached_thoraxe_input_dir::Union{Nothing, AbstractString},
        specieslist_filter::Bool,
        biomart_datasets_filter::Bool;
        specieslist_resolver::Function = _resolve_effective_specieslist,
        biomart_resolver::Function = _resolve_biomart_datasets_specieslist)
    species_filter = if cached_thoraxe_input_dir === nothing && specieslist_filter
        specieslist_resolver(target, specieslist, orthology)
    else
        (specieslist = _normalized_specieslist(specieslist), warnings = String[])
    end
    biomart_filter = if cached_thoraxe_input_dir === nothing && biomart_datasets_filter
        biomart_resolver(target, species_filter.specieslist)
    else
        (specieslist = species_filter.specieslist, warnings = String[])
    end
    return (;
        species_filter,
        biomart_filter,
        effective_specieslist = biomart_filter.specieslist
    )
end
