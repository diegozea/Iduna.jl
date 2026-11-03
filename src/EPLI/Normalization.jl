"""
    comparable_positions_normalization(score; reference_score=nothing)

Normalize validated identity counts as `100 * matched_positions / comparable_positions`.
The explicit `(0, 0)` convention returns zero. Per-sample records retain the counts
and whether any positions were compared, independently of the chosen normalizer.
"""
function comparable_positions_normalization(score; reference_score = nothing)
    counts = _normalization_counts(score)
    comparable = counts.comparable_positions
    comparable == 0 && return NormalizationResult(0.0, 0.0)
    return NormalizationResult(100 * Float64(counts.matched_positions) / comparable, comparable)
end

_normalization_counts(counts::IdentityCounts) = counts
_normalization_counts(score::AbstractScore) = _require_identity_counts(score)
function _normalization_counts(score::NamedTuple)
    IdentityCounts(score.matched_positions, score.comparable_positions)
end

"""
    self_reference_normalization(score; reference_score)

Normalize a raw score by the score from aligning the reference MSA to itself.
Both raw values must be finite; the reference score must also be nonzero.
"""
function self_reference_normalization(score; reference_score = nothing)
    reference_score === nothing &&
        error("self_reference_normalization requires a reference self score.")
    numerator = raw_score(_score_result(score))
    denominator = raw_score(_score_result(reference_score))
    denominator == 0.0 &&
        error("Cannot use self-reference normalization when the reference self score is zero.")
    return NormalizationResult(100 * numerator / denominator, denominator)
end

"""
    no_normalization(score; reference_score=nothing)

Use the finite raw score directly as the EPLI component, without a denominator.
"""
function no_normalization(score; reference_score = nothing)
    return NormalizationResult(raw_score(_score_result(score)))
end

"""
    requires_reference_score(normalization_fn)

Whether a normalizer needs a reference self score. Custom callable normalizers
can extend this method; callers may also set `compute_reference_self_score` explicitly.
"""
function requires_reference_score(normalization_fn)
    normalization_fn === self_reference_normalization
end
