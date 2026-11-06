"""
    EPLIScore(value)

A validated EPLI value. Its immutable `value::Float64` is finite; negative
values and values above 100 are allowed for alternative score functions.
Conversion to Float64 precedes validation. Invalid values raise `ArgumentError`.
An absent score is represented by `nothing` in the containing result.
"""
struct EPLIScore
    value::Float64
    function EPLIScore(value::Real)
        return new(_finite_score(value, "EPLI"))
    end
end

function _finite_score(value::Real, label::AbstractString)
    converted = Float64(value)
    # isfinite returns false for NaN, Inf, and -Inf.
    isfinite(converted) ||
        throw(ArgumentError("$(label) must be finite; got $(converted)."))
    return converted
end

"""
    epli_value(result)

Return the EPLI as a `Float64` from an `EPLI.epli_score` result or a `SeedSelection`.
Raise `ArgumentError` if the score is unavailable, for example when seed sampling
was disabled.
"""
epli_value(result) = _require_epli(result.epli)

_require_epli(score::EPLIScore) = score.value
function _require_epli(::Nothing)
    throw(ArgumentError("An EPLI score is required but unavailable."))
end

_optional_epli(value::Real) = EPLIScore(value)
function _optional_epli(value::Union{Nothing, Missing, EPLIScore})
    return ismissing(value) ? nothing : value
end

function _validate_epli_reason(score::Union{Nothing, EPLIScore},
        reason::Union{Nothing, Symbol})
    reason === nothing && return nothing
    score === nothing ||
        throw(ArgumentError("An available EPLI cannot have an absence reason."))
    reason in (:sampling_disabled, :ineligible, :not_recorded) ||
        throw(ArgumentError("Unknown EPLI absence reason: $(reason)."))
    return reason
end
