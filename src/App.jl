import ArgParse
import JSON3

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args

const _APP_TIMEOUT_KEYS = Set([
    :transcript_query_timeout_seconds,
    :transcript_query_timeout_max_seconds,
    :thoraxe_timeout_seconds
])

function _parse_pid_thresholds(value::AbstractString)
    vals = parse.(Float64, strip.(split(value, ',')))
    isempty(vals) && error("--pid-thresholds cannot be empty.")
    return vals
end

function _parse_timeout_value(value::AbstractString)
    normalized = lowercase(strip(value))
    normalized in ("none", "nothing", "off", "0") && return nothing
    parsed = parse(Float64, value)
    parsed > 0 || error("Timeout values must be positive, or use none.")
    return parsed
end

function _throw_parse_error(settings, err)
    throw(err)
end

function _app_arg_settings(;
        exc_handler::Function = _throw_parse_error,
        exit_after_help::Bool = false)
    settings = ArgParseSettings(;
        prog = "iduna",
        description = "Build and expand one ThorAxe-based MSA from a UniProt accession or an Ensembl transcript ID.",
        usage = "iduna <UniProt-or-Ensembl-transcript-ID> --mmseqs-db <path> [options]",
        autofix_names = true,
        exc_handler,
        exit_after_help,
        help_alignment_width = 36
    )
    @add_arg_table! settings begin
        "id"
        help = "UniProt accession or Ensembl transcript ID."
        metavar = "UniProt-or-Ensembl-transcript-ID"
        arg_type = String
        required = true

        "--mmseqs-db"
        help = "MMseqs2 database prefix."
        metavar = "path"
        arg_type = String
        required = true

        "--workdir"
        help = "Output/work directory."
        metavar = "path"
        arg_type = String

        "--output-dir"
        help = "Alias for --workdir."
        metavar = "path"
        arg_type = String

        "--overwrite"
        help = "Regenerate package-owned outputs in workdir."
        action = :store_true

        "--centroids"
        help = "Also save the centroid-level MSA before MMseqs2 cluster expansion."
        action = :store_true

        "--pid-thresholds"
        help = "Comma-separated thresholds, default 10,20,30,60,80."
        metavar = "list"
        arg_type = String

        "--transcript-id"
        help = "Ensembl transcript used to disambiguate UniProt input."
        metavar = "id"
        arg_type = String

        "--ensembl-gene-id"
        help = "Explicit Ensembl gene ID."
        metavar = "id"
        arg_type = String

        "--ensembl-protein-id"
        help = "Explicit Ensembl protein ID."
        metavar = "id"
        arg_type = String

        "--species"
        help = "Species name passed to transcript_query."
        metavar = "name"
        arg_type = String

        "--specieslist"
        help = "Species list passed to ThorAxe."
        metavar = "path-or-list"
        arg_type = String

        "--orthology"
        help = "Orthology relationship passed to transcript_query: 1:1, 1:n, or m:n."
        metavar = "relationship"
        arg_type = String

        "--no-specieslist-filter"
        help = "Skip Ensembl homology specieslist filtering."
        dest_name = "specieslist_filter"
        action = :store_false

        "--no-biomart-datasets-filter"
        help = "Skip BioMart dataset availability specieslist filtering."
        dest_name = "biomart_datasets_filter"
        action = :store_false

        "--thoraxe-input-dir"
        help = "Reuse a complete transcript_query bundle."
        metavar = "path"
        arg_type = String

        "--transcript-query-timeout-seconds"
        help = "Stop transcript_query after n seconds."
        metavar = "n"
        arg_type = String

        "--transcript-query-timeout-max-seconds"
        help = "Maximum timeout after retry backoff."
        metavar = "n"
        arg_type = String

        "--transcript-query-retries"
        help = "Number of transcript_query attempts."
        metavar = "n"
        arg_type = Int

        "--no-specieslist-timeout-fallback"
        help = "Keep specieslist after transcript_query timeout."
        dest_name = "allow_specieslist_timeout_fallback"
        action = :store_false

        "--thoraxe-timeout-seconds"
        help = "Stop each thoraxe run after n seconds."
        metavar = "n"
        arg_type = String

        "--threads"
        help = "Threads for MMseqs2."
        metavar = "n"
        arg_type = Int
    end
    return settings
end

function _postprocess_app_args(parsed::Dict{Symbol, Any})
    kwargs = Dict{Symbol, Any}()
    for (key, value) in parsed
        value === nothing && continue
        # ArgParse reads a few options as strings so users can pass "none".
        if key === :pid_thresholds
            kwargs[key] = _parse_pid_thresholds(value)
        elseif key in _APP_TIMEOUT_KEYS
            kwargs[key] = _parse_timeout_value(value)
        else
            kwargs[key] = value
        end
    end
    return kwargs
end

function _parse_app_args(args::Vector{String};
        exc_handler::Function = _throw_parse_error,
        exit_after_help::Bool = false)
    parsed = parse_args(args, _app_arg_settings(; exc_handler, exit_after_help);
        as_symbols = true)
    parsed === nothing && return Dict{Symbol, Any}(:help => true)
    return _postprocess_app_args(parsed)
end

function _run_app(args::Vector{String} = ARGS)
    kwargs = _parse_app_args(args; exc_handler = ArgParse.default_handler)
    if get(kwargs, :help, false)
        return nothing
    end
    result = iduna(; kwargs...)
    JSON3.pretty(stdout, Utils.result_summary(result))
    println()
    return result
end

function (@main)(args)
    _run_app(String.(args))
    return nothing
end
