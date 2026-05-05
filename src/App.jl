using JSON3

function app_help()
    return """
    Usage:
      iduna <UniProt-or-Ensembl-transcript-ID> --mmseqs-db <path> [options]

    Required:
      --mmseqs-db <path>              MMseqs2 database prefix.

    Options:
      --workdir <path>                Output/work directory.
      --output-dir <path>             Alias for --workdir.
      --overwrite                     Regenerate package-owned outputs in workdir.
      --pid-thresholds <list>         Comma-separated thresholds, default 10,20,30,60,80.
      --transcript-id <id>            Ensembl transcript used to disambiguate UniProt input.
      --ensembl-gene-id <id>          Explicit Ensembl gene ID.
      --ensembl-protein-id <id>       Explicit Ensembl protein ID.
      --species <name>                Species name passed to transcript_query.
      --specieslist <path-or-list>    Species list passed to ThorAxe.
      --thoraxe-input-dir <path>      Reuse a complete transcript_query bundle.
      --transcript-query-timeout-seconds <n>
                                      Stop transcript_query after n seconds.
      --transcript-query-timeout-max-seconds <n>
                                      Maximum timeout after retry backoff.
      --transcript-query-retries <n>  Number of transcript_query attempts.
      --no-specieslist-timeout-fallback
                                      Keep specieslist after transcript_query timeout.
      --thoraxe-timeout-seconds <n>   Stop each thoraxe run after n seconds.
      --threads <n>                   Threads for MMseqs2.
      --help                          Show this message.
    """
end

function parse_pid_thresholds(value::AbstractString)
    vals = parse.(Float64, strip.(split(value, ',')))
    isempty(vals) && error("--pid-thresholds cannot be empty.")
    return vals
end

function parse_timeout_value(value::AbstractString)
    normalized = lowercase(strip(value))
    normalized in ("none", "nothing", "off", "0") && return nothing
    parsed = parse(Float64, value)
    parsed > 0 || error("Timeout values must be positive, or use none.")
    return parsed
end

function parse_app_args(args::Vector{String})
    kwargs = Dict{Symbol, Any}()
    positional = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            kwargs[:help] = true
            i += 1
        elseif arg == "--overwrite"
            kwargs[:overwrite] = true
            i += 1
        elseif arg == "--no-specieslist-timeout-fallback"
            kwargs[:allow_specieslist_timeout_fallback] = false
            i += 1
        elseif startswith(arg, "--")
            key = Symbol(replace(arg[3:end], '-' => '_'))
            i += 1
            i <= length(args) || error("Missing value for $(arg).")
            value = args[i]
            if key === :pid_thresholds
                kwargs[key] = parse_pid_thresholds(value)
            elseif key === :threads || key === :transcript_query_retries
                kwargs[key] = parse(Int, value)
            elseif key === :transcript_query_timeout_seconds ||
                   key === :transcript_query_timeout_max_seconds ||
                   key === :thoraxe_timeout_seconds
                kwargs[key] = parse_timeout_value(value)
            else
                kwargs[key] = value
            end
            i += 1
        else
            push!(positional, arg)
            i += 1
        end
    end
    if get(kwargs, :help, false)
        return kwargs
    end
    length(positional) <= 1 || error("Pass only one positional input ID.")
    !isempty(positional) && (kwargs[:id] = first(positional))
    haskey(kwargs, :mmseqs_db) || error("--mmseqs-db is required.")
    return kwargs
end

function run_app(args::Vector{String} = ARGS)
    kwargs = parse_app_args(args)
    if get(kwargs, :help, false)
        print(app_help())
        return nothing
    end
    result = iduna(; kwargs...)
    JSON3.pretty(stdout, Utils.result_summary(result))
    println()
    return result
end

function (@main)(args)
    run_app(String.(args))
    return nothing
end
