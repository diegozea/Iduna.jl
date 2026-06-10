"""
    comparable_positions_normalization(score; reference_score=nothing)

Normalize an HHsuite identity score as `100 * matched_positions / comparable_positions`.
"""
function comparable_positions_normalization(score; reference_score = nothing)
    comparable = Int(score.comparable_positions)
    comparable == 0 && return (; normalized_score = 0.0, normalization_score = 0.0)
    return (;
        normalized_score = 100 * Float64(score.matched_positions) / comparable,
        normalization_score = Float64(comparable))
end

"""
    self_reference_normalization(score; reference_score)

Normalize a score by the raw score obtained from aligning the reference MSA to itself.
"""
function self_reference_normalization(score; reference_score = nothing)
    reference_score === nothing &&
        error("self_reference_normalization requires a reference self score.")
    denominator = Float64(reference_score.raw_score)
    denominator == 0.0 &&
        error("Cannot use self-reference normalization when the reference self score is zero.")
    return (;
        normalized_score = 100 * Float64(score.raw_score) / denominator,
        normalization_score = denominator)
end

"""
    no_normalization(score; reference_score=nothing)

Use `score.raw_score` directly as the EPLI component.
"""
function no_normalization(score; reference_score = nothing)
    return (;
        normalized_score = Float64(score.raw_score),
        normalization_score = missing)
end

function _normalization_result(value)
    if value isa Real
        return (; normalized_score = Float64(value), normalization_score = missing)
    elseif value isa NamedTuple
        :normalized_score in keys(value) ||
            error("Normalization functions that return a named tuple must include `normalized_score`.")
        normalization_score = :normalization_score in keys(value) ?
                              value.normalization_score : missing
        return (;
            normalized_score = Float64(value.normalized_score),
            normalization_score)
    end
    error("Normalization function must return a real number or a named tuple.")
end
