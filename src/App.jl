import ArgParse
import JSON

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args

function _parse_pid_thresholds(value::AbstractString)
    vals = parse.(Float64, strip.(split(value, ',')))
    isempty(vals) && error("--pid-thresholds cannot be empty.")
    return vals
end

function _throw_parse_error(settings, err)
    throw(err)
end

function _app_arg_settings(;
        exc_handler::Function = _throw_parse_error,
        exit_after_help::Bool = false)
    settings = ArgParseSettings(;
        prog = "iduna",
        description = "Build one ThorAxe-based MSA from a UniProt accession or an Ensembl transcript ID, " *
                      "optionally expanding it with MMseqs2/HMMER.",
        usage = "iduna <UniProt-or-Ensembl-transcript-ID> [--mmseqs-db <path>] [--no-expansion] [options]",
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
        help = "MMseqs2 database prefix. Required unless --no-expansion is used."
        metavar = "path"
        arg_type = String

        "--no-expansion"
        help = "Stop after the ThorAxe MSA stage without MMseqs2/HMMER expansion."
        action = :store_true

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

        "--transcript-query-retries"
        help = "Number of transcript_query attempts."
        metavar = "n"
        arg_type = Int

        "--pid-sample-count"
        help = "PID-specific species samples per candidate MSA for seed selection, default 45."
        metavar = "n"
        arg_type = Int

        "--pid-sample-fraction"
        help = "Fraction of non-reference species retained in each PID-specific sample, default 0.8."
        metavar = "x"
        arg_type = Float64

        "--pid-sample-seed"
        help = "Random seed for PID species sampling. A random seed is recorded if omitted."
        metavar = "n"
        arg_type = Int

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
        if key === :pid_thresholds
            kwargs[key] = _parse_pid_thresholds(value)
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

function _run_app(args::Vector{String} = ARGS; runner::Function = iduna)
    kwargs = _parse_app_args(args; exc_handler = ArgParse.default_handler)
    if get(kwargs, :help, false)
        return nothing
    end
    result = runner(; kwargs...)
    JSON.json(stdout, Utils.result_summary(result); pretty = true)
    println()
    return result
end

function (@main)(args)
    _run_app(String.(args))
    return nothing
end
