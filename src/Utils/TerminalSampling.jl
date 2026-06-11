_io_is_tty(io) = io isa Base.TTY
_io_is_tty(io::IOContext) = _io_is_tty(getfield(io, :io))

function _env_truthy(name::AbstractString)
    value = lowercase(strip(get(ENV, name, "")))
    return value in ("1", "true", "yes", "on")
end

function _terminal_progress_enabled(output = stderr)
    (_env_truthy("CI") || _env_truthy("GITHUB_ACTIONS")) && return false
    return _io_is_tty(output)
end

function _sample_rng(seed::UInt64, sample_idx::Integer)
    # Give each sample its own repeatable random stream from the same run seed.
    mixed = xor(seed, UInt64(sample_idx) * 0xbf58476d1ce4e5b9)
    return MersenneTwister(Int(mod(mixed, UInt64(typemax(Int)))))
end

function _sample_indices(n_total::Integer, reference_idx::Integer,
        fraction::Real, rng::MersenneTwister)
    # Always keep the reference sequence so every sample remains comparable.
    selectable = [i for i in 1:n_total if i != reference_idx]
    isempty(selectable) && return [reference_idx]
    n_keep = clamp(
        round(Int, Float64(fraction) * length(selectable)), 1, length(selectable))
    return vcat(reference_idx, sample(rng, selectable, n_keep; replace = false))
end
