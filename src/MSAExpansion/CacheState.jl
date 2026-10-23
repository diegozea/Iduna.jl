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

function _validate_expansion_alignment(path::AbstractString, format)
    if format === Stockholm
        text = strip(read(path, String))
        startswith(text, r"# STOCKHOLM 1\.0\r?\n") && endswith(text, "//") ||
            error("Expansion Stockholm file is incomplete: $(path).")
    end
    alignment = read_file(path, format; keepinserts = true)
    nsequences(alignment) > 0 && ncolumns(alignment) > 0 ||
        error("Expansion alignment has no sequences or columns: $(path).")
    return nothing
end

function _validate_expansion_outputs(outputs::NamedTuple)
    for (name, path) in _required_expansion_outputs(outputs)
        if endswith(name, "hits_fasta")
            # No hits is a valid result, but a present record must contain residues.
            filesize(path) == 0 && continue
            records = read_file(path, FASTASequences)
            !isempty(records) && all(sequence -> length(sequence) > 0, records) ||
                error("Expansion hit FASTA has no complete sequences: $(path).")
        else
            format = endswith(name, "stockholm") ? Stockholm : A3M
            _validate_expansion_alignment(path, format)
        end
    end
    return nothing
end

function _write_step_state(run_dir::AbstractString,
        status::Symbol,
        identity,
        outputs::NamedTuple;
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing,
        action::Union{Nothing, Symbol, AbstractString} = nothing)
    status === :done && _validate_expansion_outputs(outputs)
    return _write_stage_state(run_dir;
        stage = "msa_expansion",
        stage_key = _expansion_stage_key(identity),
        status,
        identity,
        outputs,
        required_outputs = _required_expansion_outputs(outputs),
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
    state = _classify_stage_state(run_dir, identity, _required_expansion_outputs(outputs);
        stage_label = "MSA expansion")
    if state.reusable && !_expansion_reuse_trusted(identity)
        return (; reusable = false, status = :stale,
            warning = "MSA expansion provenance is incomplete; rebuilding cached outputs.")
    end
    return state
end

function _write_cache_warning(logs_dir::AbstractString, warning::AbstractString)
    mkpath(logs_dir)
    open(joinpath(logs_dir, "cache_warning.log"), "a") do io
        println(io, "[", now(UTC), "] ", warning)
    end
    return nothing
end
