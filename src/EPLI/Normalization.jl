"""
    comparable_positions_normalization(score; reference_score=nothing)

Normalize validated identity counts as `100 * matched_positions / comparable_positions`.
The explicit `(0, 0)` convention returns zero. Per-sample records retain the counts
and whether any positions were compared, independently of the chosen normalizer.
"""
function comparable_positions_normalization(score::ScoreResult;
        reference_score::Union{Nothing, ScoreResult} = nothing)
    counts = _require_identity_counts(score)
    comparable = counts.comparable_positions
    comparable == 0 && return NormalizationResult(0.0, 0.0)
    return NormalizationResult(100 * Float64(counts.matched_positions) / comparable, comparable)
end

"""
    self_reference_normalization(score; reference_score)

Normalize a raw score by the score from aligning the reference MSA to itself.
Both raw values must be finite; the reference score must also be nonzero.
"""
function self_reference_normalization(score::ScoreResult;
        reference_score::Union{Nothing, ScoreResult} = nothing)
    reference_score === nothing &&
        error("self_reference_normalization requires a reference self score.")
    numerator = raw_score(score)
    denominator = raw_score(reference_score)
    denominator == 0.0 &&
        error("Cannot use self-reference normalization when the reference self score is zero.")
    return NormalizationResult(100 * numerator / denominator, denominator)
end

"""
    no_normalization(score; reference_score=nothing)

Use the finite raw score directly as the EPLI component, without a denominator.
"""
function no_normalization(score::ScoreResult;
        reference_score::Union{Nothing, ScoreResult} = nothing)
    return NormalizationResult(raw_score(score))
end

"""
    requires_reference_score(normalization_fn)

Whether a normalizer needs a reference self score. Custom callable normalizers
can extend this method; callers may also set `compute_reference_self_score` explicitly.
"""
function requires_reference_score(normalization_fn)
    normalization_fn === self_reference_normalization
end
