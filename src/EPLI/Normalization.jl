function _validate_finite_score(value::Float64, label::AbstractString)
    # isfinite returns false for NaN, Inf, and -Inf.
    isfinite(value) || throw(ArgumentError("$(label) must be finite; got $(value)."))
    return nothing
end

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

Normalize an HHsuite identity score as `100 * matched_positions / comparable_positions`.
"""
function comparable_positions_normalization(score; reference_score = nothing)
    comparable = Int(score.comparable_positions)
    matched = Float64(score.matched_positions)
    _validate_finite_score(matched, "Matched positions")
    comparable == 0 && return (; normalized_score = 0.0, normalization_score = 0.0)
    return _normalization_result((;
        normalized_score = _percentage_ratio(matched, Float64(comparable)),
        normalization_score = Float64(comparable)))
end

"""
    self_reference_normalization(score; reference_score)

Normalize a score by the raw score obtained from aligning the reference MSA to itself.
"""
function self_reference_normalization(score; reference_score = nothing)
    reference_score === nothing &&
        error("self_reference_normalization requires a reference self score.")
    raw_score = Float64(score.raw_score)
    denominator = Float64(reference_score.raw_score)
    _validate_finite_score(raw_score, "Raw score")
    _validate_finite_score(denominator, "Reference self score")
    denominator == 0.0 &&
        error("Cannot use self-reference normalization when the reference self score is zero.")
    return _normalization_result((;
        normalized_score = _percentage_ratio(raw_score, denominator),
        normalization_score = denominator))
end

"""
    no_normalization(score; reference_score=nothing)

Use `score.raw_score` directly as the EPLI component.
"""
function no_normalization(score; reference_score = nothing)
    raw_score = Float64(score.raw_score)
    _validate_finite_score(raw_score, "Raw score")
    return (;
        normalized_score = raw_score,
        normalization_score = missing)
end

function _normalization_result(value)
    if value isa Real
        normalized_score = Float64(value)
        normalization_score = missing
    elseif value isa NamedTuple
        :normalized_score in keys(value) ||
            error("Normalization functions that return a named tuple must include `normalized_score`.")
        normalization_score = :normalization_score in keys(value) ?
                              value.normalization_score : missing
        normalized_score = Float64(value.normalized_score)
        normalization_score = ismissing(normalization_score) ? missing :
                              Float64(normalization_score)
    else
        error("Normalization function must return a real number or a named tuple.")
    end
    _validate_finite_score(normalized_score, "Normalized score")
    ismissing(normalization_score) ||
        _validate_finite_score(normalization_score, "Normalization score")
    return (; normalized_score, normalization_score)
end
