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

abstract type AbstractScore end

"""
    IdentityScore(counts::IdentityCounts; diagnostics=(;))

A profile comparison measured by identity counts. `raw_score(score)` is derived
from the matched count, so the raw score and counts cannot disagree. Optional
tool diagnostics are retained separately and are not score invariants.
"""
struct IdentityScore{D} <: AbstractScore
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
struct RawScore{D} <: AbstractScore
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

Return the finite raw value of an `IdentityScore` or `RawScore`.
"""
raw_score(score::IdentityScore) = Float64(score.counts.matched_positions)
raw_score(score::RawScore) = score.raw_score

function _require_identity_counts(score::AbstractScore)
    score.counts === nothing && throw(ArgumentError("This score has no identity counts."))
    return score.counts
end

"""
    NormalizationResult(value, normalization_score=missing;
                        no_comparable_positions=false)

An immutable normalized component with finite Float64 values. A missing
`normalization_score` means no denominator was supplied, as with `no_normalization`.
Negative values and ratios above 100 are valid. `no_comparable_positions=true`
records the explicit zero-comparison convention and requires both values to be zero.
Invalid values raise `ArgumentError` at construction, not on field access.
"""
struct NormalizationResult
    normalized_score::Float64
    normalization_score::Union{Missing, Float64}
    no_comparable_positions::Bool
    function NormalizationResult(value::Real, denominator::Union{Missing, Real} = missing;
            no_comparable_positions::Bool = false)
        normalized = _finite_score(value, "Normalized score")
        scale = ismissing(denominator) ? missing :
                _finite_score(denominator, "Normalization score")
        if no_comparable_positions
            if ismissing(scale) || !iszero(scale) || !iszero(normalized)
                throw(ArgumentError("No comparable positions requires a zero component and denominator."))
            end
        end
        return new(normalized, scale, no_comparable_positions)
    end
end

function _callback_identity_counts(value::NamedTuple)
    matched = get(value, :matched_positions, missing)
    comparable = get(value, :comparable_positions, missing)
    (ismissing(matched) || ismissing(comparable)) && return nothing
    return IdentityCounts(matched, comparable)
end

_score_result(value::Union{RawScore, IdentityScore}) = value
function _score_result(value::NamedTuple)
    return RawScore(value.raw_score;
        counts = _callback_identity_counts(value), diagnostics = value)
end

function _normalization_result(value)
    value isa NormalizationResult && return value
    value isa Real && return NormalizationResult(value)
    value isa NamedTuple ||
        error("Normalization function must return NormalizationResult, a real number, or a named tuple.")
    :normalized_score in keys(value) ||
        error("Normalization functions that return a named tuple must include `normalized_score`.")
    return NormalizationResult(
        value.normalized_score, get(value, :normalization_score, missing);
        no_comparable_positions = get(value, :no_comparable_positions, false))
end
