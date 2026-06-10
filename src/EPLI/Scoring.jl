_function_label(fn) = string(fn)

function _score_field(score, name::Symbol)
    name in keys(score) || return missing
    return getfield(score, name)
end

function _score_sample(spec,
        reference_msa_fasta::AbstractString,
        reference_score;
        score_fn::Function,
        normalization_fn::Function,
        logs_dir::Union{Nothing, AbstractString},
        overwrite::Bool)
    build_sample_msa = spec.build_sample_msa
    sample_msa_fasta = build_sample_msa(spec; overwrite)
    isfile(sample_msa_fasta) ||
        error("Sample MSA FASTA is missing: $(sample_msa_fasta).")
    score = score_fn(reference_msa_fasta, sample_msa_fasta;
        logs_dir = logs_dir === nothing ? nothing : joinpath(logs_dir, "hhalign"),
        label = spec.sample_label)
    normalized = _normalization_result(normalization_fn(score; reference_score))
    n_reference = spec.n_sequences_reference
    n_sample = spec.n_sequences_sample === nothing ?
               nsequences(read_file(sample_msa_fasta, FASTA)) :
               spec.n_sequences_sample
    return (;
        sample = spec.sample,
        sample_label = spec.sample_label,
        raw_score = Float64(score.raw_score),
        normalization_score = normalized.normalization_score,
        normalized_score = normalized.normalized_score,
        epli_component = normalized.normalized_score,
        n_sequences_reference = n_reference,
        n_sequences_sample = n_sample,
        reference_msa_fasta,
        sample_msa_fasta,
        matched_positions = _score_field(score, :matched_positions),
        comparable_positions = _score_field(score, :comparable_positions),
        hhalign_score = _score_field(score, :hhalign_score),
        probability = _score_field(score, :probability),
        e_value = _score_field(score, :e_value))
end

function _reference_self_score(reference_msa_fasta::AbstractString,
        compute_reference_self_score::Bool;
        score_fn::Function,
        logs_dir::Union{Nothing, AbstractString})
    compute_reference_self_score || return nothing
    return score_fn(reference_msa_fasta, reference_msa_fasta;
        logs_dir = logs_dir === nothing ? nothing : joinpath(logs_dir, "hhalign"),
        label = "reference_self")
end

function _score_sample_specs!(rows::AbstractVector,
        sample_specs::AbstractVector,
        reference_msa_fasta::AbstractString,
        reference_score;
        score_fn::Function,
        normalization_fn::Function,
        logs_dir::Union{Nothing, AbstractString},
        overwrite::Bool,
        parallel::Bool,
        progress_desc::AbstractString,
        progress_output::IO,
        progress_enabled::Bool)
    run_in_parallel = parallel && length(sample_specs) > 1 && Threads.threadpoolsize() > 1
    progress = ProgressMeter.Progress(length(sample_specs);
        desc = progress_desc,
        output = progress_output,
        enabled = progress_enabled)
    ProgressMeter.update!(progress, 0; force = true)
    progress_lock = ReentrantLock()
    function advance_progress!()
        lock(progress_lock) do
            ProgressMeter.next!(progress)
        end
        return nothing
    end
    if run_in_parallel
        Threads.@threads :greedy for idx in eachindex(sample_specs)
            rows[idx] = _score_sample(sample_specs[idx], reference_msa_fasta,
                reference_score; score_fn, normalization_fn, logs_dir, overwrite)
            advance_progress!()
        end
    else
        for idx in eachindex(sample_specs)
            rows[idx] = _score_sample(sample_specs[idx], reference_msa_fasta,
                reference_score; score_fn, normalization_fn, logs_dir, overwrite)
            advance_progress!()
        end
    end
    return rows
end

function _score_alignment_samples(reference_msa_fasta::AbstractString,
        sample_specs::AbstractVector;
        score_fn::Function = hhsuite_identity_score,
        normalization_fn::Function = comparable_positions_normalization,
        scores_path::Union{Nothing, AbstractString} = nothing,
        logs_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        compute_reference_self_score::Bool = normalization_fn ===
                                             self_reference_normalization,
        progress_desc::AbstractString = "Scoring EPLI samples: ",
        progress_output::IO = stderr,
        progress_enabled::Bool = _terminal_progress_enabled(progress_output),
        parallel::Bool = true)
    isempty(sample_specs) && error("At least one sample is required to compute EPLI.")
    reference_score = _reference_self_score(
        reference_msa_fasta, compute_reference_self_score;
        score_fn, logs_dir)
    rows = Vector{Union{Nothing, NamedTuple}}(nothing, length(sample_specs))
    _score_sample_specs!(rows, sample_specs, reference_msa_fasta, reference_score;
        score_fn,
        normalization_fn,
        logs_dir,
        overwrite,
        parallel,
        progress_desc,
        progress_output,
        progress_enabled)
    score_rows = [row for row in rows if row !== nothing]
    isempty(score_rows) && error("No sample scores were computed.")
    if scores_path !== nothing
        mkpath(dirname(scores_path))
        CSV.write(scores_path, DataFrame(score_rows))
    end
    return (;
        rows = score_rows,
        epli = median([row.epli_component for row in score_rows]),
        reference_score,
        scores_path)
end

"""
    epli_score(input_fasta, workdir, aligner_fn; aligner_args=Cmd(String[]))

Compute EPLI for an aligner over sequence-row subsamples from an unaligned FASTA.

# Arguments

- `input_fasta::AbstractString`: unaligned FASTA with the full sequence set.
- `workdir::AbstractString`: output directory for inputs, MSAs, logs, and tables.
- `aligner_fn::Function`: callable that writes an aligned FASTA from an input FASTA.
  It must accept `logs_dir`, `run_label`, and `aligner_args` keywords.

# Keywords

- `score_fn = hhsuite_identity_score`: pairwise profile score function.
- `normalization_fn = comparable_positions_normalization`: score normalization.
- `aligner_args::Cmd = Cmd(String[])`: extra command arguments passed to `aligner_fn`.
- `sample_count::Integer = 45`: number of sequence-row samples.
- `sample_fraction::Real = 0.8`: fraction of non-reference rows per sample.
- `sample_seed = nothing`: random seed; `nothing` records a generated seed.
- `reference_sequence = nothing`: reference sequence name, or first row when `nothing`.
- `overwrite::Bool = false`: rebuild existing output files.
"""
function epli_score(input_fasta::AbstractString,
        workdir::AbstractString,
        aligner_fn::Function;
        score_fn::Function = hhsuite_identity_score,
        normalization_fn::Function = comparable_positions_normalization,
        aligner_args::Cmd = Cmd(String[]),
        sample_count::Integer = 45,
        sample_fraction::Real = 0.8,
        sample_seed::Union{Nothing, Integer} = nothing,
        reference_sequence = nothing,
        overwrite::Bool = false,
        compute_reference_self_score::Bool = normalization_fn ===
                                             self_reference_normalization,
        progress_output::IO = stderr,
        progress_enabled::Bool = _terminal_progress_enabled(progress_output))
    _validate_sampling_options(sample_count, sample_fraction)
    seed = sample_seed === nothing ? UInt64(rand(UInt32)) :
           _normalize_sample_seed(sample_seed)
    fasta = _fasta_records(input_fasta)
    reference_idx = _reference_index(fasta.names, reference_sequence)
    full_input = joinpath(workdir, "input", "full_sequences.fasta")
    reference_msa = joinpath(workdir, "msa", "reference_msa.fasta")
    logs_dir = joinpath(workdir, "logs")
    scores_path = joinpath(workdir, "scores.csv")
    summary_path = joinpath(workdir, "summary.csv")
    write_fasta(full_input, fasta.records)
    reference_msa = _run_aligner(aligner_fn, full_input, reference_msa;
        logs_dir = joinpath(logs_dir, "aligner"),
        run_label = "full",
        overwrite,
        aligner_args)
    sample_specs = NamedTuple[]
    for sample_idx in 1:Int(sample_count)
        label = _sample_label(sample_idx)
        rng = _sample_rng(seed, sample_idx)
        indices = _sample_indices(length(fasta.records), reference_idx, sample_fraction, rng)
        sample_input = joinpath(workdir, "samples", "$(label)_sequences.fasta")
        sample_msa = joinpath(workdir, "samples", "$(label)_msa.fasta")
        write_fasta(sample_input, [fasta.records[i] for i in indices])
        push!(sample_specs,
            (;
                sample = sample_idx,
                sample_label = label,
                sample_sequence_fasta = sample_input,
                sample_msa_fasta = sample_msa,
                n_sequences_reference = length(fasta.records),
                n_sequences_sample = length(indices),
                build_sample_msa = (spec;
                    overwrite = false) -> _run_aligner(
                    aligner_fn, spec.sample_sequence_fasta, spec.sample_msa_fasta;
                    logs_dir = joinpath(logs_dir, "aligner"),
                    run_label = spec.sample_label,
                    overwrite,
                    aligner_args)))
    end
    scored = _score_alignment_samples(reference_msa, sample_specs;
        score_fn,
        normalization_fn,
        scores_path,
        logs_dir,
        overwrite,
        compute_reference_self_score,
        progress_desc = "Scoring EPLI samples: ",
        progress_output,
        progress_enabled)
    summary = (;
        epli = scored.epli,
        n_samples = length(scored.rows),
        sample_count = Int(sample_count),
        sample_fraction = Float64(sample_fraction),
        sample_seed = seed,
        score_function = _function_label(score_fn),
        normalization_function = _function_label(normalization_fn),
        aligner_args = string(aligner_args),
        input_fasta = String(input_fasta),
        reference_msa_fasta = reference_msa,
        scores_path)
    CSV.write(summary_path, DataFrame([summary]))
    return merge(summary,
        (;
            workdir = String(workdir),
            summary_path,
            rows = scored.rows,
            reference_score = scored.reference_score))
end
