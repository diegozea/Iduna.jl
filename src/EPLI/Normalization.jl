function _percentage_ratio(numerator::Float64, denominator::Float64)
    scaled = 100 * numerator
    value = scaled / denominator
    if isfinite(scaled) && floatmin(Float64) <= abs(value) < floatmax(Float64) / 2
        return value
    end
    # Avoid intermediate overflow and double rounding near the Float64 limits.
    # Julia 1.12 scopes these settings to the current task, including parallel scoring.
    return setprecision(BigFloat, 256) do
        setrounding(BigFloat, RoundNearest) do
            Float64(100 * BigFloat(numerator) / BigFloat(denominator), RoundNearest)
        end
    end
end

"""
    comparable_positions_normalization(score; reference_score=nothing)

Normalize validated identity counts as `100 * matched_positions / comparable_positions`.
The explicit `(0, 0)` convention returns zero and records `no_comparable_positions=true`.
"""
function comparable_positions_normalization(score; reference_score = nothing)
    counts = _normalization_counts(score)
    comparable = counts.comparable_positions
    comparable == 0 && return NormalizationResult(0.0, 0.0; no_comparable_positions = true)
    return NormalizationResult(
        _percentage_ratio(Float64(counts.matched_positions), Float64(comparable)),
        comparable)
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
    return NormalizationResult(_percentage_ratio(numerator, denominator), denominator)
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
requires_reference_score(::Any) = false
requires_reference_score(::typeof(self_reference_normalization)) = true
