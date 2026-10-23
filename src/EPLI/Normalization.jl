function _finite_score(value, label::AbstractString)
    score = Float64(value)
    isfinite(score) || throw(ArgumentError("$(label) must be finite; got $(value)."))
    return score
end

"""
    comparable_positions_normalization(score; reference_score=nothing)

Normalize an HHsuite identity score as `100 * matched_positions / comparable_positions`.
"""
function comparable_positions_normalization(score; reference_score = nothing)
    comparable = Int(score.comparable_positions)
    matched = _finite_score(score.matched_positions, "Matched positions")
    comparable == 0 && return (; normalized_score = 0.0, normalization_score = 0.0)
    return (;
        normalized_score = _finite_score(100 * matched / comparable, "Normalized score"),
        normalization_score = Float64(comparable))
end

"""
    self_reference_normalization(score; reference_score)

Normalize a score by the raw score obtained from aligning the reference MSA to itself.
"""
function self_reference_normalization(score; reference_score = nothing)
    reference_score === nothing &&
        error("self_reference_normalization requires a reference self score.")
    denominator = _finite_score(reference_score.raw_score, "Reference self score")
    denominator == 0.0 &&
        error("Cannot use self-reference normalization when the reference self score is zero.")
    return (;
        normalized_score = _finite_score(
            100 * _finite_score(score.raw_score, "Raw score") / denominator,
            "Normalized score"),
        normalization_score = denominator)
end

"""
    no_normalization(score; reference_score=nothing)

Use `score.raw_score` directly as the EPLI component.
"""
function no_normalization(score; reference_score = nothing)
    return (;
        normalized_score = _finite_score(score.raw_score, "Raw score"),
        normalization_score = missing)
end

function _normalization_result(value)
    if value isa Real
        return (; normalized_score = _finite_score(value, "Normalized score"),
            normalization_score = missing)
    elseif value isa NamedTuple
        :normalized_score in keys(value) ||
            error("Normalization functions that return a named tuple must include `normalized_score`.")
        normalization_score = :normalization_score in keys(value) ?
                              value.normalization_score : missing
        normalization_score = ismissing(normalization_score) ? missing :
                              _finite_score(normalization_score, "Normalization score")
        return (;
            normalized_score = _finite_score(value.normalized_score, "Normalized score"),
            normalization_score)
    end
    error("Normalization function must return a real number or a named tuple.")
end
