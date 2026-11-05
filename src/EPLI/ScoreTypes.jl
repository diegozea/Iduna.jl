"""
    IdentityCounts(matched_positions, comparable_positions)

Immutable integer counts satisfying `0 ≤ matched_positions ≤ comparable_positions`.
Integer-valued real inputs are accepted. Impossible, non-finite, fractional, or
out-of-range counts raise `ArgumentError`. `(0, 0)` is valid and normalizes to zero.
"""
struct IdentityCounts
    matched_positions::Int
    comparable_positions::Int
    function IdentityCounts(matched::Real, comparable::Real)
        0 <= matched <= comparable ||
            throw(ArgumentError("Identity counts require 0 ≤ matched ≤ comparable."))
        return new(_identity_count(matched), _identity_count(comparable))
    end
end

function _identity_count(value::Real)
    isfinite(value) && isinteger(value) && 0 <= value <= typemax(Int) ||
        throw(ArgumentError("Identity counts must be nonnegative Int values; got $(value)."))
    return Int(value)
end

"""
    ScoreResult

Abstract type for results returned by a score callback. Subtypes must implement
`raw_score`, returning a real value representable as a finite Float64. They may
also implement `identity_counts` and `score_diagnostics`.

Accessors must return consistent values during scoring and be safe for concurrent
reads. Iduna checks raw values and counts when callbacks return, before normalization.
"""
abstract type ScoreResult end

"""
    IdentityScore(counts::IdentityCounts; diagnostics=(;))

A profile comparison measured by identity counts. `raw_score(score)` is derived
from the matched count, so the raw score and counts cannot disagree. Optional
tool diagnostics are retained separately and are not score invariants.
"""
struct IdentityScore{D} <: ScoreResult
    counts::IdentityCounts
    diagnostics::D
end
function IdentityScore(counts::IdentityCounts; diagnostics = (;))
    IdentityScore(counts, diagnostics)
end

"""
    RawScore(value; counts=nothing, diagnostics=(;))

A finite Float64 raw score, including negative profile or custom scores.
Optional `counts` must be validated `IdentityCounts`. Diagnostics may contain
missing or non-finite values when they are not used as the scientific score.
"""
struct RawScore{D} <: ScoreResult
    raw_score::Float64
    counts::Union{Nothing, IdentityCounts}
    diagnostics::D
    function RawScore(value::Real; counts::Union{Nothing, IdentityCounts} = nothing,
            diagnostics = (;))
        return new{typeof(diagnostics)}(_finite_score(value, "Raw score"), counts, diagnostics)
    end
end

"""
    raw_score(score)

Return the raw numerical value of a score. Every `ScoreResult` subtype must
implement this method. The value must be real and representable as a finite Float64.
For `IdentityScore`, it is the matched count; for `RawScore`, it is the stored value.
"""
raw_score(score::IdentityScore) = Float64(score.counts.matched_positions)
raw_score(score::RawScore) = score.raw_score

"""
    identity_counts(score::ScoreResult)

Return validated `IdentityCounts`, or `nothing` when the score has no identity
counts. The default is `nothing`. Implement this method to support
`comparable_positions_normalization` and identity-count columns in sample output.
"""
identity_counts(::ScoreResult) = nothing
identity_counts(score::Union{IdentityScore, RawScore}) = score.counts

"""
    score_diagnostics(score::ScoreResult)

Return optional information from the scoring tool. The default is an empty named
tuple. Prefer named fields such as `(; e_value=1e-6)`; `hhalign_score`, `probability`,
and `e_value` properties are copied to sample output. Other information remains
available to custom normalizers. These fields do not supply identity counts or
change the raw score, and may contain missing or non-finite values.
"""
score_diagnostics(::ScoreResult) = (;)
score_diagnostics(score::Union{IdentityScore, RawScore}) = score.diagnostics

# Check extension results once at the callback boundary, retaining values for output.
function _validated_score_values(score::ScoreResult)
    raw = _finite_score(raw_score(score)::Real, "Raw score")
    counts = identity_counts(score)::Union{Nothing, IdentityCounts}
    return raw, counts
end

function _require_identity_counts(score::ScoreResult)
    counts = identity_counts(score)::Union{Nothing, IdentityCounts}
    counts === nothing && throw(ArgumentError("This score has no identity counts."))
    return counts
end

"""
    NormalizationResult(value, normalization_score=missing)

An immutable normalized component with finite Float64 values. A missing
`normalization_score` means no denominator was supplied, as with `no_normalization`.
Negative values and ratios above 100 are valid. Identity counts belong to the
original score; normalization does not determine whether positions were compared.
Invalid values raise `ArgumentError` at construction, not on field access.
"""
struct NormalizationResult
    normalized_score::Float64
    normalization_score::Union{Missing, Float64}
    function NormalizationResult(value::Real, denominator::Union{Missing, Real} = missing)
        normalized = _finite_score(value, "Normalized score")
        scale = ismissing(denominator) ? missing :
                _finite_score(denominator, "Normalization score")
        return new(normalized, scale)
    end
end
