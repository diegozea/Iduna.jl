# Cache State
# -----------

function _exception_summary(err)
    return (;
        type = string(typeof(err)),
        message = sprint(showerror, err)
    )
end

function _expansion_stage_key(identity)
    return "expansion:$(identity.target.ensembl_gene_id):$(identity.target.transcript_id):$(format_pid_dir(identity.seed.pid))"
end

function _write_step_state(run_dir::AbstractString,
        status::Symbol,
        identity,
        outputs::NamedTuple;
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing,
        action::Union{Nothing, Symbol, AbstractString} = nothing)
    return _write_stage_state(run_dir;
        stage = "msa_expansion",
        stage_key = _expansion_stage_key(identity),
        status,
        identity,
        outputs,
        warnings,
        exception,
        action,
        extra = (; step = "msa_expansion"))
end

_read_step_state(run_dir::AbstractString) = _read_stage_state(run_dir)

function _step_state_unreadable_message(state)
    if state === nothing
        return "state file disappeared while reading"
    elseif state isa NamedTuple && haskey(state, :unreadable)
        return state.unreadable
    end
    return nothing
end

function _classify_step_state(run_dir::AbstractString, identity, outputs::NamedTuple)
    return _classify_stage_state(run_dir, identity, _required_expansion_outputs(outputs);
        stage_label = "MSA expansion")
end

function _write_cache_warning(logs_dir::AbstractString, warning::AbstractString)
    mkpath(logs_dir)
    open(joinpath(logs_dir, "cache_warning.log"), "a") do io
        println(io, "[", now(UTC), "] ", warning)
    end
    return nothing
end
