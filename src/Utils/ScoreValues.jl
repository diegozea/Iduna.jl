"""
    ComputedEPLI(value)

An EPLI value that was computed. Its immutable `value::Float64` is finite;
negative values and values above 100 are allowed for alternative score functions.
Conversion to Float64 precedes validation. Invalid values raise `ArgumentError`.
"""
struct ComputedEPLI
    value::Float64
    function ComputedEPLI(value::Real)
        return new(_finite_score(value, "EPLI"))
    end
end

"""
    UncomputedEPLI(reason)

An EPLI that was not computed. `reason` is `:sampling_disabled`, `:ineligible`,
or `:not_recorded` for an older result without score-state information.
This state cannot be used as a numeric score or ranked.
"""
struct UncomputedEPLI
    reason::Symbol
    function UncomputedEPLI(reason::Symbol)
        reason in (:sampling_disabled, :ineligible, :not_recorded) ||
            throw(ArgumentError("Unknown uncomputed EPLI reason: $(reason)."))
        return new(reason)
    end
end

function _finite_score(value::Real, label::AbstractString)
    converted = Float64(value)
    # isfinite returns false for NaN, Inf, and -Inf.
    isfinite(converted) ||
        throw(ArgumentError("$(label) must be finite; got $(converted)."))
    return converted
end

_epli_value(score::ComputedEPLI) = score.value
_epli_value(::UncomputedEPLI) = missing
function _epli_state(score::Union{ComputedEPLI, UncomputedEPLI})
    score isa ComputedEPLI ? "computed" : "not_computed"
end
_epli_reason(::ComputedEPLI) = nothing
_epli_reason(score::UncomputedEPLI) = String(score.reason)

"""
    require_epli(score)

Return a computed EPLI value, or throw `ArgumentError` if it was not computed.
"""
require_epli(score::ComputedEPLI) = score.value
function require_epli(score::UncomputedEPLI)
    throw(ArgumentError("EPLI was not computed: $(score.reason)."))
end

_epli_score(value::Real) = ComputedEPLI(value)
_epli_score(::Missing) = UncomputedEPLI(:not_recorded)
_epli_score(value::Union{ComputedEPLI, UncomputedEPLI}) = value
