using Base.Threads
import JSON

"""
    iduna(; id, mmseqs_db=nothing, no_expansion=false, kwargs...) -> IdunaResult
    iduna(id; mmseqs_db=nothing, no_expansion=false, kwargs...) -> IdunaResult

Build one ThorAxe-based MSA from a UniProt accession or an Ensembl transcript
ID and, by default, expand it with MMseqs2/HMMER. Set `no_expansion=true` to
stop after the ThorAxe MSA stage without requiring an MMseqs2 database.

The function is file-first: every external tool writes its inputs, outputs, and
logs under `workdir`, and the returned [`IdunaResult`](@ref) stores paths plus
the resolved identifiers and validation statistics.

`orthology` controls the ThorAxe relationship filter (`"1:1"`, `"1:n"`, or
`"m:n"`). By default, Iduna starts from `specieslist="ases"`, the Ases
webserver's 12-species list. Set `specieslist="all"` or `specieslist=""` for
unrestricted ThorAxe species selection, or pass a comma-separated species list,
file path, or single species name. The lowercase strings `"ases"` and `"all"`
are reserved presets; use `./ases` or `./all` for files with those names. Iduna
filters the requested species with Ensembl homology and then filters against
currently available BioMart Ensembl Gene datasets before `transcript_query`.
That step can be slow, especially for broad species lists; runs are often faster
when a small curated `specieslist` is provided. `pid_sample_count` controls how
many PID-specific species samples are drawn from each candidate `msa_0` and
scored against that PID's full-species MSA; the default is 45.
`pid_sample_fraction` controls the fraction of non-reference species retained in
each sample, and `pid_sample_seed` can make the sampling reproducible. If
omitted, a random seed is recorded with the result. `sampling_strategy` controls
how species samples are shared across PID thresholds; the default `:common`
uses one shared set of species samples for all eligible PID candidates.
Candidate `msa_0`
reconstructions with indels versus UniProt are reported and excluded from seed
selection; substitution-only differences are warnings. Seed selection uses
highest median identity, then highest mean identity, then the largest candidate
`msa_0` species count, then the first PID in `pid_thresholds` order. Set
`pid_sample_count=0` to skip seed selection and carry every eligible PID
candidate forward. Pass `thoraxe_input_dir` to reuse a complete
transcript_query bundle instead of fetching it again. Set `centroids=true` to
also save the centroid-level MSA before MMseqs2 expands centroid hits to cluster
members; the main expansion and validation still use the full expanded MSA.
`centroids=true` requires expansion and cannot be combined with
`no_expansion=true`. To run the repeated ThorAxe sample runs faster, start Julia
with more than one Julia thread, for example `julia --threads 4`; Iduna uses
those threads automatically during PID sample scoring. The `mmseqs_threads`
keyword controls MMseqs2 expansion only.
"""
function iduna(; id::Union{Nothing, AbstractString} = nothing,
        uniprot_id::Union{Nothing, AbstractString} = nothing,
        ensembl_transcript_id::Union{Nothing, AbstractString} = nothing,
        transcript_id::Union{Nothing, AbstractString} = nothing,
        ensembl_gene_id::Union{Nothing, AbstractString} = nothing,
        ensembl_protein_id::Union{Nothing, AbstractString} = nothing,
        mmseqs_db::Union{Nothing, AbstractString} = nothing,
        no_expansion::Bool = false,
        workdir::Union{Nothing, AbstractString} = nothing,
        output_dir::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        pid_thresholds::AbstractVector{<:Real} = Utils.DEFAULT_PID_THRESHOLDS,
        species::Union{Nothing, AbstractString} = nothing,
        specieslist::AbstractString = "ases",
        orthology::AbstractString = "1:1",
        specieslist_filter::Bool = true,
        biomart_datasets_filter::Bool = true,
        thoraxe_input_dir::Union{Nothing, AbstractString} = nothing,
        transcript_query_retries::Integer = 2,
        pid_sample_count::Integer = 45,
        pid_sample_fraction::Real = 0.8,
        pid_sample_seed::Union{Nothing, Integer} = nothing,
        sampling_strategy::Symbol = :common,
        match_mode::Integer = 1,
        match_ratio::Union{Nothing, Real} = nothing,
        hmmbuild_symfrac::Real = 0.0,
        centroids::Bool = false,
        mmseqs_threads::Union{Nothing, Integer} = Threads.nthreads(),
        _resolve_target::Function = IDMapping.resolve_target,
        _build_thoraxe_msa::Function = ThorAxeMSA.build_thoraxe_msa,
        _expand_msa::Function = MSAExpansion.expand_msa,
        _validate_results::Function = ResultsValidation.validate_results)
    primary, disambiguating_transcript,
    supplied_uniprot = _normalize_primary_input(;
        id, uniprot_id, ensembl_transcript_id, transcript_id)
    _validate_expansion_options(; mmseqs_db, no_expansion, centroids)
    root = nothing
    target = nothing
    thoraxe = nothing
    expansions = Union{Nothing, Utils.ExpansionResult}[]
    validations = Utils.ValidationResult[]
    failed_stage = "prepare_output_dir"

    try
        root = Utils.prepare_output_dir(primary; workdir, output_dir, overwrite)
        prepared_details = overwrite ?
                           (; input_id = primary, workdir = root, overwrite) :
                           (; input_id = primary, workdir = root)
        @info "Prepared Iduna output directory." prepared_details...

        failed_stage = "resolve_target"
        target = _resolve_target_stage(primary, root;
            supplied_uniprot,
            ensembl_gene_id,
            ensembl_protein_id,
            disambiguating_transcript,
            species,
            overwrite,
            resolver = _resolve_target)

        failed_stage = "thoraxe_msa"
        @info "Starting ThorAxe pipeline."
        thoraxe = _build_thoraxe_msa(target, root;
            pid_thresholds,
            specieslist,
            orthology,
            specieslist_filter,
            biomart_datasets_filter,
            cached_thoraxe_input_dir = thoraxe_input_dir,
            overwrite,
            transcript_query_retries,
            pid_sample_count,
            pid_sample_fraction,
            pid_sample_seed,
            sampling_strategy)

        if !no_expansion
            failed_stage = "msa_expansion"
            n_seeds = length(thoraxe.seeds)
            for (index, seed) in enumerate(thoraxe.seeds)
                @info "Expanding MSA seed." seed_index=index n_seeds
                push!(expansions,
                    _expand_msa(target, seed, root;
                        mmseqs_db,
                        overwrite,
                        match_mode,
                        match_ratio,
                        hmmbuild_symfrac,
                        centroids,
                        mmseqs_threads))
            end
        else
            @info "Skipping MSA expansion."
        end

        failed_stage = "validation"
        @info "Validating Iduna results." n_seeds=length(thoraxe.seeds)
        for (index, seed) in enumerate(thoraxe.seeds)
            push!(validations,
                _call_validate_results(_validate_results, target, seed,
                    no_expansion ? nothing : expansions[index], root; overwrite))
        end
        warnings = vcat(target.warnings, thoraxe.warnings,
            Iterators.flatten(validation.warnings for validation in validations)...)
        status = _pipeline_status(warnings)
        current_stage_keys = _current_stage_keys(target, thoraxe, expansions, validations)
        current_stages = Utils.collect_stage_summaries(root;
            stage_keys = current_stage_keys)
        result_identity = (; input_id = primary,
            stage_hashes = [stage.identity_hash
                            for stage in current_stages
                            if stage.identity_hash !== nothing &&
                               stage.stage_key != "result"])
        Utils._write_stage_state(Utils._pipeline_stage_dir(root, "result");
            stage = "result",
            stage_key = "result",
            status = :done,
            identity = result_identity,
            outputs = (; result_json = joinpath(root, "result.json")),
            action = :run,
            workdir = root)
        stages = Utils.collect_stage_summaries(root;
            stage_keys = _current_stage_keys(target, thoraxe, expansions, validations;
                include_result = true))
        result = Utils.IdunaResult(;
            input_id = primary,
            workdir = root,
            target,
            thoraxe_msa = thoraxe,
            expansions,
            validations,
            stages,
            warnings,
            status
        )
        result = Utils._relative_result_paths(result)
        @info "Writing Iduna result artifact." result_path=joinpath(
            root, "result.json")
        Utils.write_json(joinpath(root, "result.json"), Utils.result_summary(result))
        @info "Iduna pipeline completed." status=result.status
        return result
    catch err
        if root !== nothing
            _write_failure_result(
                primary, root, failed_stage, err; target, thoraxe, expansions,
                validations)
        end
        rethrow()
    end
end

iduna(id::AbstractString; kwargs...) = iduna(; id, kwargs...)

function _validate_expansion_options(; mmseqs_db, no_expansion::Bool, centroids::Bool)
    if no_expansion && centroids
        error("centroids=true requires MMseqs expansion; use no_expansion=false.")
    end
    if !no_expansion && mmseqs_db === nothing
        error("Pass mmseqs_db for expansion, or set no_expansion=true to stop after the ThorAxe MSA stage.")
    end
    return nothing
end

function _pipeline_status(warnings::AbstractVector{<:AbstractString})
    isempty(warnings) ? :ok : :warn
end

function _call_validate_results(validator::Function,
        target,
        seed,
        expansion,
        workdir;
        overwrite::Bool)
    if validator === ResultsValidation.validate_results
        return validator(target, seed, expansion, workdir; overwrite)
    end
    return validator(target, seed, expansion, workdir)
end

function _pid_stage_key(prefix::AbstractString, pid::Real)
    return "$(prefix):$(Utils.format_pid_dir(pid))"
end

function _expansion_stage_key(target, pid::Real)
    return "expansion:$(target.ensembl_gene_id):$(target.transcript_id):$(Utils.format_pid_dir(pid))"
end

function _current_stage_count(count::Integer, total::Integer,
        failed_stage::Union{Nothing, AbstractString},
        stage_name::AbstractString)
    return failed_stage == stage_name && count < total ? count + 1 : count
end

function _push_expansion_stage_keys!(keys::Vector{String},
        target,
        seeds,
        count::Integer)
    target === nothing && return keys
    for seed in Iterators.take(seeds, min(count, length(seeds)))
        push!(keys, _expansion_stage_key(target, seed.pid))
    end
    return keys
end

function _push_validation_stage_keys!(keys::Vector{String},
        seeds,
        count::Integer)
    for seed in Iterators.take(seeds, min(count, length(seeds)))
        push!(keys, _pid_stage_key("validation", seed.pid))
    end
    return keys
end

function _current_stage_keys(target,
        thoraxe,
        expansions::AbstractVector,
        validations::AbstractVector;
        include_result::Bool = false,
        failed_stage::Union{Nothing, AbstractString} = nothing)
    keys = String["target"]
    if target !== nothing || thoraxe !== nothing
        push!(keys, "thoraxe_input", "thoraxe_msa")
    end
    if thoraxe !== nothing
        total = length(thoraxe.seeds)
        n_expansion_keys = _current_stage_count(
            length(expansions), total, failed_stage, "msa_expansion")
        n_validation_keys = _current_stage_count(
            length(validations), total, failed_stage, "validation")
        _push_expansion_stage_keys!(keys, target, thoraxe.seeds, n_expansion_keys)
        _push_validation_stage_keys!(keys, thoraxe.seeds, n_validation_keys)
    end
    include_result && push!(keys, "result")
    return unique(keys)
end

function _write_failure_result(input_id::AbstractString, workdir::AbstractString,
        failed_stage::AbstractString, err; target = nothing, thoraxe = nothing,
        expansions = Union{Nothing, Utils.ExpansionResult}[],
        validations = Utils.ValidationResult[])
    try
        Utils.write_json(joinpath(workdir, "result.json"),
            _failure_result_summary(input_id, workdir, failed_stage, err;
                target, thoraxe, expansions, validations))
    catch write_err
        @warn "Could not write failure result artifact." workdir exception=(write_err,
            catch_backtrace())
    end
    return nothing
end

function _failure_expansions(failed_stage::AbstractString, thoraxe,
        expansions::AbstractVector)
    padded = Union{Nothing, Utils.ExpansionResult}[expansions...]
    if String(failed_stage) == "msa_expansion" && thoraxe !== nothing
        append!(padded, fill(nothing, max(0, length(thoraxe.seeds) - length(padded))))
    end
    return padded
end

function _failure_result_summary(input_id::AbstractString, workdir::AbstractString,
        failed_stage::AbstractString, err; target = nothing, thoraxe = nothing,
        expansions = Union{Nothing, Utils.ExpansionResult}[],
        validations = Utils.ValidationResult[])
    stage_keys = _current_stage_keys(target, thoraxe, expansions, validations;
        failed_stage)
    expansions = _failure_expansions(failed_stage, thoraxe, expansions)
    warnings = _partial_warnings(target, thoraxe, validations)
    if target !== nothing && thoraxe !== nothing
        partial = Utils.IdunaResult(;
            input_id = String(input_id),
            workdir = String(workdir),
            target,
            thoraxe_msa = thoraxe,
            expansions,
            validations,
            stages = Utils.collect_stage_summaries(workdir; stage_keys),
            warnings,
            status = :error)
        summary = Utils.result_summary(partial)
        return merge(summary,
            (;
                failed_stage = String(failed_stage),
                exception = _exception_summary(err, workdir)))
    end
    return (;
        input_id = String(input_id),
        status = "error",
        failed_stage = String(failed_stage),
        stages = Utils.collect_stage_summaries(workdir; stage_keys),
        warnings,
        target = target === nothing ? nothing : _target_summary(target, workdir),
        exception = _exception_summary(err, workdir)
    )
end

function _partial_warnings(target, thoraxe, validations)
    warnings = String[]
    target !== nothing && append!(warnings, target.warnings)
    thoraxe !== nothing && append!(warnings, thoraxe.warnings)
    if validations !== nothing
        for validation in validations
            append!(warnings, validation.warnings)
        end
    end
    return warnings
end

function _exception_summary(err, _workdir::Union{Nothing, AbstractString} = nothing)
    return (;
        type = string(typeof(err)),
        message = sprint(showerror, err)
    )
end

function _provided_primary_ids(; id, uniprot_id, ensembl_transcript_id)
    provided_primary = Pair{Symbol, String}[]
    id !== nothing && push!(provided_primary, :id => String(id))
    ensembl_transcript_id !== nothing &&
        push!(provided_primary, :ensembl_transcript_id => String(ensembl_transcript_id))
    if isempty(provided_primary) && uniprot_id !== nothing
        push!(provided_primary, :uniprot_id => String(uniprot_id))
    end
    return provided_primary
end

function _transcript_primary_input(transcript_id)
    transcript_id === nothing &&
        error("Pass id, uniprot_id, ensembl_transcript_id, or transcript_id.")
    return String(transcript_id), nothing, nothing
end

function _consistent_primary_value(provided_primary)
    first_value = first(provided_primary).second
    if any(pair -> pair.second != first_value, provided_primary)
        error("Conflicting primary identifiers were provided: $(provided_primary).")
    end
    return first_value
end

function _validate_uniprot_alias(kind::Symbol, first_value::String, uniprot_id)
    if kind === :uniprot && uniprot_id !== nothing && String(uniprot_id) != first_value
        error("Conflicting UniProt identifiers were provided: id=$(first_value), uniprot_id=$(uniprot_id).")
    end
    return nothing
end

function _normalize_primary_input(; id, uniprot_id, ensembl_transcript_id, transcript_id)
    # Accept old and new argument names, but still require one clear input.
    provided_primary = _provided_primary_ids(; id, uniprot_id, ensembl_transcript_id)
    isempty(provided_primary) && return _transcript_primary_input(transcript_id)

    first_value = _consistent_primary_value(provided_primary)
    kind = Utils.id_kind(first_value)
    # For UniProt input, transcript_id is only a tie-breaker among candidates.
    _validate_uniprot_alias(kind, first_value, uniprot_id)
    disambiguating_transcript = kind === :uniprot ? transcript_id : nothing
    supplied_uniprot = kind === :uniprot ? first_value :
                       (uniprot_id === nothing ? nothing : String(uniprot_id))
    return first_value, disambiguating_transcript, supplied_uniprot
end

function _target_identity(primary::AbstractString;
        supplied_uniprot,
        ensembl_gene_id,
        ensembl_protein_id,
        disambiguating_transcript,
        species)
    return (;
        input_id = String(primary),
        input_kind = String(Utils.id_kind(primary)),
        uniprot_id = supplied_uniprot === nothing ? nothing : String(supplied_uniprot),
        ensembl_gene_id = ensembl_gene_id === nothing ? nothing : String(ensembl_gene_id),
        ensembl_protein_id = ensembl_protein_id === nothing ? nothing :
                             String(ensembl_protein_id),
        transcript_id = disambiguating_transcript === nothing ? nothing :
                        String(disambiguating_transcript),
        species = species === nothing ? nothing : String(species)
    )
end

_target_stage_dir(workdir::AbstractString) = Utils._pipeline_stage_dir(workdir, "target")

function _target_from_summary(data, workdir::AbstractString)
    return Utils.ResolvedTarget(;
        input_id = String(data["input_id"]),
        input_kind = Symbol(String(data["input_kind"])),
        uniprot_id = get(data, "uniprot_id", nothing),
        ensembl_gene_id = String(data["ensembl_gene_id"]),
        transcript_id = String(data["transcript_id"]),
        ensembl_protein_id = get(data, "ensembl_protein_id", nothing),
        species = get(data, "species", nothing),
        uniprot_sequence_path = get(data, "uniprot_sequence_path", nothing),
        ensembl_protein_sequence_path = get(data, "ensembl_protein_sequence_path", nothing),
        sequence_validated = get(data, "sequence_validated", nothing),
        mapping_confirmed = get(data, "mapping_confirmed", nothing),
        workdir,
        warnings = String.(get(data, "warnings", String[]))
    )
end

function _target_artifacts_exist(target::Utils.ResolvedTarget, workdir::AbstractString)
    for path in (target.uniprot_sequence_path, target.ensembl_protein_sequence_path)
        path === nothing && continue
        isfile(Utils._resolve_artifact_path(path, workdir)) || return false
    end
    return true
end

function _read_cached_target(workdir::AbstractString)
    path = joinpath(workdir, "target.json")
    isfile(path) || return nothing
    try
        target = _target_from_summary(JSON.parse(read(path, String)), workdir)
        _target_artifacts_exist(target, workdir) || return nothing
        return target
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _clear_target_outputs!(workdir::AbstractString)
    target_json = joinpath(workdir, "target.json")
    isfile(target_json) && rm(target_json; force = true)
    sequences_dir = joinpath(workdir, "sequences")
    isdir(sequences_dir) && Utils.safe_rm(sequences_dir, workdir)
    return nothing
end

function _resolve_target_stage(primary::AbstractString, workdir::AbstractString;
        supplied_uniprot,
        ensembl_gene_id,
        ensembl_protein_id,
        disambiguating_transcript,
        species,
        overwrite::Bool,
        resolver::Function)
    identity = _target_identity(primary;
        supplied_uniprot,
        ensembl_gene_id,
        ensembl_protein_id,
        disambiguating_transcript,
        species)
    outputs = (; target_json = joinpath(workdir, "target.json"))
    stage_dir = _target_stage_dir(workdir)
    cache = overwrite ? (; reusable = false, status = :stale, warning = nothing) :
            Utils._classify_stage_state(stage_dir, identity, outputs;
        stage_label = "target")
    cached_target = cache.reusable ? _read_cached_target(workdir) : nothing
    if cached_target !== nothing
        @info "Reusing cached target metadata."
        Utils._write_stage_state(stage_dir;
            stage = "target",
            stage_key = "target",
            status = :done,
            identity,
            outputs,
            action = :reuse,
            workdir,
            warnings = cached_target.warnings)
        return cached_target
    end

    cache.warning === nothing || @warn String(cache.warning) workdir status=cache.status
    action = cache.status === :missing ? :run : :rebuild
    (overwrite || cache.status !== :missing) && _clear_target_outputs!(workdir)
    Utils._write_stage_state(stage_dir;
        stage = "target",
        stage_key = "target",
        status = :running,
        identity,
        outputs,
        action,
        workdir)
    try
        @info "Resolving target identifiers."
        target = resolver(primary;
            workdir,
            uniprot_id = supplied_uniprot,
            ensembl_gene_id,
            ensembl_protein_id,
            transcript_id = disambiguating_transcript,
            species)
        @info "Writing target metadata." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id
        Utils.write_json(outputs.target_json, _target_summary(target, workdir))
        Utils._write_stage_state(stage_dir;
            stage = "target",
            stage_key = "target",
            status = :done,
            identity,
            outputs,
            action,
            workdir,
            warnings = target.warnings)
        return target
    catch err
        Utils._write_stage_state(stage_dir;
            stage = "target",
            stage_key = "target",
            status = :failed,
            identity,
            outputs,
            action,
            workdir,
            exception = _exception_summary(err, workdir))
        rethrow()
    end
end

_summary_path(path, workdir::Nothing) = path
_summary_path(path, workdir::AbstractString) = Utils._relative_artifact_path(path, workdir)

function _target_summary(target::Utils.ResolvedTarget,
        workdir::Union{Nothing, AbstractString} = nothing)
    return (;
        input_id = target.input_id,
        input_kind = String(target.input_kind),
        uniprot_id = target.uniprot_id,
        ensembl_gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        ensembl_protein_id = target.ensembl_protein_id,
        species = target.species,
        uniprot_sequence_path = _summary_path(target.uniprot_sequence_path, workdir),
        ensembl_protein_sequence_path = _summary_path(
            target.ensembl_protein_sequence_path, workdir),
        sequence_validated = target.sequence_validated,
        mapping_confirmed = target.mapping_confirmed,
        warnings = target.warnings
    )
end

"""
    load_result(workdir) -> IdunaResult

Reconstruct an [`IdunaResult`](@ref) from an existing Iduna output directory
without running pipeline stages or writing files. The loader expects the current
move-safe result schema, where `result.json` contains structured ThorAxe seed
records and artifact paths are relative to `workdir`.
"""
function load_result(workdir::AbstractString)
    root = abspath(String(workdir))
    summary = _load_result_summary(root)
    status = _result_status(summary)
    status === :error &&
        @warn "Loading partial IdunaResult from a failed run." workdir=root failed_stage=get(
            summary, "failed_stage", nothing)
    target = _load_result_target(summary, root)
    target === nothing &&
        error("Cannot load IdunaResult from $(root): target metadata is unavailable.")
    thoraxe = _load_result_thoraxe(summary, target, root, status)
    expansions = _load_result_expansions(summary, target, thoraxe, root)
    validations = _load_result_validations(summary, thoraxe, expansions, root)
    result = Utils.IdunaResult(;
        input_id = String(get(summary, "input_id", target.input_id)),
        workdir = root,
        target,
        thoraxe_msa = thoraxe,
        expansions,
        validations,
        stages = collect(Any, _result_vector(summary, "stages")),
        warnings = String.(get(summary, "warnings", String[])),
        status
    )
    return Utils._relative_result_paths(result)
end

function _load_result_summary(workdir::AbstractString)
    path = joinpath(workdir, "result.json")
    isfile(path) ||
        error("Cannot load IdunaResult from $(workdir): result.json is missing.")
    data = JSON.parse(read(path, String))
    data isa AbstractDict ||
        error("Cannot load IdunaResult from $(workdir): result.json is not an object.")
    return data
end

_result_status(summary) = Symbol(String(get(summary, "status", "ok")))

function _result_section(summary, key::AbstractString)
    value = get(summary, key, nothing)
    return value isa AbstractDict ? value : nothing
end

function _result_vector(summary, key::AbstractString)
    get(summary, key, Any[]) isa AbstractVector ? get(summary, key, Any[]) : Any[]
end

_path_or_default(::Nothing, default::AbstractString, _workdir::AbstractString) = default
function _path_or_default(path, _default::AbstractString, workdir::AbstractString)
    return Utils._resolve_artifact_path(String(path), workdir)
end

function _nullable_path_or_default(data::AbstractDict, key::AbstractString,
        default::AbstractString, workdir::AbstractString)
    haskey(data, key) || return default
    value = data[key]
    value === nothing && return nothing
    return _path_or_default(value, default, workdir)
end

_maybe_int_result(value) = value === nothing ? nothing : Int(value)
_maybe_float_result(value) = value === nothing ? nothing : Float64(value)
_maybe_bool_result(value) = value === nothing ? nothing : Bool(value)
_maybe_string_result(value) = value === nothing ? nothing : String(value)
_missing_or_float_result(value) = value === nothing ? missing : Float64(value)

function _target_from_result_summary(summary, data, workdir::AbstractString)
    input_id = String(get(data, "input_id", get(summary, "input_id", "")))
    input_kind = haskey(data, "input_kind") ? Symbol(String(data["input_kind"])) :
                 Utils.id_kind(input_id)
    return Utils.ResolvedTarget(;
        input_id,
        input_kind,
        uniprot_id = get(data, "uniprot_id", nothing),
        ensembl_gene_id = String(data["ensembl_gene_id"]),
        transcript_id = String(data["transcript_id"]),
        ensembl_protein_id = get(data, "ensembl_protein_id", nothing),
        species = get(data, "species", nothing),
        uniprot_sequence_path = get(data, "uniprot_sequence_path", nothing),
        ensembl_protein_sequence_path = get(data, "ensembl_protein_sequence_path", nothing),
        sequence_validated = get(data, "sequence_validated", nothing),
        mapping_confirmed = get(data, "mapping_confirmed", nothing),
        workdir,
        warnings = String.(get(data, "warnings", String[]))
    )
end

function _load_result_target(summary, workdir::AbstractString)
    target_json = joinpath(workdir, "target.json")
    if isfile(target_json)
        try
            return _target_from_summary(JSON.parse(read(target_json, String)), workdir)
        catch err
            err isa InterruptException && rethrow()
            @warn "Could not read target.json; falling back to result.json target metadata." workdir target_json exception=(
                err, catch_backtrace())
        end
    end
    data = _result_section(summary, "target")
    data === nothing && return nothing
    @warn "Loading target metadata from result.json; sequence artifact paths may be incomplete." workdir
    try
        return _target_from_result_summary(summary, data, workdir)
    catch err
        err isa InterruptException && rethrow()
        @warn "Could not reconstruct target metadata from result.json." workdir exception=(
            err, catch_backtrace())
        return nothing
    end
end

function _load_result_seed_record(data, workdir::AbstractString)
    data isa AbstractDict || error("ThorAxe seed record is not an object.")
    return Utils.SeedSelection(;
        pid = Float64(data["pid"]),
        median_identity = _missing_or_float_result(get(data, "median_identity", nothing)),
        mean_identity = _missing_or_float_result(get(data, "mean_identity", nothing)),
        stockholm_path = String(data["stockholm_path"]),
        fasta_path = get(data, "fasta_path", nothing) === nothing ? nothing :
                     String(data["fasta_path"]),
        s_exon_blocks_tsv = get(data, "s_exon_blocks_tsv", nothing) === nothing ?
                            nothing : String(data["s_exon_blocks_tsv"]),
        summary_path = String(data["summary_path"]),
        used_fallback_dir = Bool(get(data, "used_fallback_dir", false)),
        workdir = String(workdir))
end

function _load_result_seeds(data, _summary_path::AbstractString, workdir::AbstractString)
    isempty(data) && return Utils.SeedSelection[]
    seeds_data = get(data, "seeds", nothing)
    seeds_data isa AbstractVector ||
        error("Cannot load IdunaResult from $(workdir): result.json uses an unsupported ThorAxe seed schema.")
    seeds = [_load_result_seed_record(seed_data, workdir) for seed_data in seeds_data]
    isempty(seeds) &&
        @warn "Loaded IdunaResult has no selected ThorAxe seeds." workdir
    return seeds
end

function _load_result_thoraxe(summary, target::Utils.ResolvedTarget,
        workdir::AbstractString, status::Symbol)
    data = _result_section(summary, "thoraxe_msa")
    if data === nothing
        @warn "result.json has no ThorAxe MSA summary; returning a partial ThorAxe result." workdir
        data = Dict{String, Any}()
    end
    summary_path = _path_or_default(get(data, "pid_summary", nothing),
        joinpath(workdir, "thoraxe_msa", "candidate_summary.csv"), workdir)
    seeds = _load_result_seeds(data, summary_path, workdir)
    artifacts = [ThorAxeMSA._seed_artifacts(seed, workdir) for seed in seeds]
    thoraxe_status = haskey(data, "status") ? Symbol(String(data["status"])) :
                     (isempty(seeds) && status === :error ? :error : :ok)
    return Utils.ThorAxeMSAResult(;
        input_dir = joinpath(workdir, "thoraxe_input"),
        thoraxe_dirs = [artifact.thoraxe_dir for artifact in artifacts],
        msa_dir = joinpath(workdir, "thoraxe_msa"),
        baseline_fastas = haskey(data, "baseline_fastas") ?
                          [Utils._resolve_artifact_path(String(path), workdir)
                           for path in data["baseline_fastas"]] :
                          [String(artifact.seed_fasta) for artifact in artifacts],
        baseline_stockholms = haskey(data, "baseline_stockholms") ?
                              [Utils._resolve_artifact_path(String(path), workdir)
                               for path in data["baseline_stockholms"]] :
                              [artifact.seed_path for artifact in artifacts],
        sequence_fastas = [artifact.sequence_fasta for artifact in artifacts],
        species_files = [artifact.species_file for artifact in artifacts],
        pid_summary = summary_path,
        seeds,
        logs_dir = joinpath(workdir, "logs", "thoraxe"),
        pid_sample_count = Int(get(data, "pid_sample_count", 0)),
        pid_sample_fraction = Float64(get(data, "pid_sample_fraction", 1.0)),
        pid_sample_seed = UInt64(get(data, "pid_sample_seed", 0)),
        sampling_strategy = Symbol(String(get(data, "sampling_strategy", "independent"))),
        warnings = String.(get(data, "warnings", String[])),
        status = thoraxe_status)
end

function _expansion_run_dir_from_summary(data, target::Utils.ResolvedTarget,
        seed::Utils.SeedSelection, workdir::AbstractString)
    run_path = get(data, "run_dir", nothing)
    run_path === nothing || return Utils._resolve_artifact_path(String(run_path), workdir)
    match_path = get(data, "match_stockholm", nothing)
    if match_path !== nothing
        expanded_dir = dirname(Utils._resolve_artifact_path(String(match_path), workdir))
        return dirname(expanded_dir)
    end
    return joinpath(workdir, "expansion", target.ensembl_gene_id,
        target.transcript_id, Utils.format_pid_dir(seed.pid))
end

function _load_result_expansion(data, target::Utils.ResolvedTarget,
        seed::Utils.SeedSelection, workdir::AbstractString)
    run_dir = _expansion_run_dir_from_summary(data, target, seed, workdir)
    outputs = MSAExpansion._expansion_output_paths(run_dir, target.transcript_id)
    seed_label = "seed_pid$(Utils.format_pid(seed.pid))"
    hits_fasta = _path_or_default(get(data, "hits_fasta", nothing),
        outputs.hits_fasta, workdir)
    seed_stockholm = Utils._resolve_artifact_path(seed.stockholm_path,
        seed.workdir === nothing ? workdir : seed.workdir)
    counts = if haskey(data, "n_hits") || !isfile(hits_fasta) || !isfile(seed_stockholm)
        (;
            n_hits = Int(get(data, "n_hits", 0)),
            n_new_hits = Int(get(data, "n_new_hits", 0)))
    else
        MSAExpansion._cached_hit_counts(hits_fasta,
            MSAExpansion._seed_id_set(seed_stockholm))
    end
    return Utils.ExpansionResult(;
        run_dir,
        seed_stockholm = _path_or_default(get(data, "seed_stockholm", nothing),
            joinpath(run_dir, "seeds", "$(seed_label).sto"), workdir),
        seed_fasta = _nullable_path_or_default(data, "seed_fasta",
            joinpath(run_dir, "seeds", "$(seed_label).fasta"), workdir),
        hits_fasta,
        full_stockholm = _path_or_default(get(data, "full_stockholm", nothing),
            outputs.full_stockholm, workdir),
        match_stockholm = _path_or_default(get(data, "match_stockholm", nothing),
            outputs.match_stockholm, workdir),
        a3m_path = _path_or_default(get(data, "a3m_path", nothing),
            outputs.a3m_path, workdir),
        s_exon_blocks_tsv = _nullable_path_or_default(data, "s_exon_blocks_tsv",
            outputs.s_exon_blocks_tsv, workdir),
        db_dir = _path_or_default(get(data, "db_dir", nothing),
            joinpath(run_dir, "dbs"), workdir),
        hmm_dir = _path_or_default(get(data, "hmm_dir", nothing),
            joinpath(run_dir, "hmm"), workdir),
        logs_dir = _path_or_default(get(data, "logs_dir", nothing),
            joinpath(run_dir, "logs"), workdir),
        n_hits = counts.n_hits,
        n_new_hits = counts.n_new_hits,
        status = Symbol(String(get(data, "status", "ok"))),
        workdir)
end

function _load_result_expansions(summary, target::Utils.ResolvedTarget,
        thoraxe::Utils.ThorAxeMSAResult, workdir::AbstractString)
    expansion_data = _result_vector(summary, "expansions")
    isempty(expansion_data) &&
        return _empty_result_expansions(summary, thoraxe, workdir)
    expansions = Union{Nothing, Utils.ExpansionResult}[_load_result_expansion_slot(
                                                           expansion_data, target,
                                                           seed, workdir, index)
                                                       for (index, seed) in
                                                           pairs(thoraxe.seeds)]
    _warn_extra_expansion_summaries(expansion_data, thoraxe.seeds, workdir)
    return expansions
end

function _empty_result_expansions(summary, thoraxe::Utils.ThorAxeMSAResult,
        workdir::AbstractString)
    failed_stage = get(summary, "failed_stage", nothing)
    expansions = Union{Nothing, Utils.ExpansionResult}[]
    failed_stage === nothing && return expansions
    String(failed_stage) == "msa_expansion" || return expansions
    isempty(thoraxe.seeds) && return expansions
    @warn "Expansion failed before producing summaries; preserving empty expansion slots." workdir n_seeds=length(
        thoraxe.seeds)
    append!(expansions, fill(nothing, length(thoraxe.seeds)))
    return expansions
end

function _load_result_expansion_slot(expansion_data::AbstractVector,
        target::Utils.ResolvedTarget, seed::Utils.SeedSelection,
        workdir::AbstractString, index::Integer)::Union{Nothing, Utils.ExpansionResult}
    data = index <= length(expansion_data) ? expansion_data[index] : nothing
    data === nothing && return _missing_result_expansion_slot(workdir, index)
    data isa AbstractDict ||
        return _malformed_result_expansion_slot(workdir, index)
    try
        return _load_result_expansion(data, target, seed, workdir)
    catch err
        err isa InterruptException && rethrow()
        @warn "Could not reconstruct expansion result; preserving empty expansion slot." workdir index exception=(
            err, catch_backtrace())
        return nothing
    end
end

function _missing_result_expansion_slot(workdir::AbstractString, index::Integer)
    @warn "Expansion summary is missing for a selected ThorAxe seed." workdir index
    return nothing
end

function _malformed_result_expansion_slot(workdir::AbstractString, index::Integer)
    @warn "Skipping partial expansion summary that cannot be matched to a ThorAxe seed." workdir index
    return nothing
end

function _warn_extra_expansion_summaries(expansion_data::AbstractVector,
        seeds::AbstractVector{Utils.SeedSelection}, workdir::AbstractString)
    length(expansion_data) <= length(seeds) && return nothing
    @warn "Ignoring expansion summaries that do not match selected ThorAxe seeds." workdir n_extra=(
        length(expansion_data) - length(seeds))
    return nothing
end

function _validation_outputs_from_summary(data, seed::Utils.SeedSelection,
        workdir::AbstractString)
    default_stats = joinpath(workdir, "validation", Utils.format_pid_dir(seed.pid),
        "stats.csv")
    stats_path = _path_or_default(get(data, "stats_path", nothing), default_stats, workdir)
    query_path = joinpath(dirname(stats_path), "query_vs_uniprot_alignment.txt")
    return (;
        stats_path,
        query_vs_uniprot_path = isfile(query_path) ? query_path : nothing)
end

function _validation_from_result_summary(data, outputs)
    return Utils.ValidationResult(;
        stats_path = outputs.stats_path,
        query_name = _maybe_string_result(get(data, "query_name", nothing)),
        query_vs_uniprot_path = outputs.query_vs_uniprot_path,
        seed_nseq = _maybe_int_result(get(data, "seed_nseq", nothing)),
        seed_ncol = _maybe_int_result(get(data, "seed_ncol", nothing)),
        seed_clusters62 = _maybe_int_result(get(data, "seed_clusters62", nothing)),
        seed_neff80 = _maybe_float_result(get(data, "seed_neff80", nothing)),
        expanded_nseq = _maybe_int_result(get(data, "expanded_nseq", nothing)),
        expanded_ncol = _maybe_int_result(get(data, "expanded_ncol", nothing)),
        expanded_clusters62 = _maybe_int_result(
            get(data, "expanded_clusters62", nothing)),
        expanded_neff80 = _maybe_float_result(get(data, "expanded_neff80", nothing)),
        aln_identical = _maybe_bool_result(get(data, "aln_identical", nothing)),
        aln_mismatches = _maybe_int_result(get(data, "aln_mismatches", nothing)),
        aln_insertions = _maybe_int_result(get(data, "aln_insertions", nothing)),
        aln_deletions = _maybe_int_result(get(data, "aln_deletions", nothing)),
        warnings = String.(get(data, "warnings", String[])),
        status = Symbol(String(get(data, "status", "ok"))))
end

function _load_result_validations(summary, thoraxe::Utils.ThorAxeMSAResult,
        expansions::Vector{Union{Nothing, Utils.ExpansionResult}},
        workdir::AbstractString)
    validation_data = _result_vector(summary, "validations")
    validations = Utils.ValidationResult[]
    for (index, data) in enumerate(validation_data)
        if !(data isa AbstractDict) || index > length(thoraxe.seeds)
            @warn "Skipping partial validation summary that cannot be matched to a ThorAxe seed." workdir index
            continue
        end
        seed = thoraxe.seeds[index]
        expansion = index <= length(expansions) ? expansions[index] : nothing
        outputs = _validation_outputs_from_summary(data, seed, workdir)
        if isfile(outputs.stats_path)
            try
                push!(validations,
                    ResultsValidation._cached_validation_result(
                        outputs, workdir, seed, expansion))
                continue
            catch err
                err isa InterruptException && rethrow()
                @warn "Could not read validation stats; falling back to result.json validation summary." workdir stats_path=outputs.stats_path exception=(
                    err, catch_backtrace())
            end
        else
            @warn "Validation stats are missing; using partial result.json validation summary." workdir stats_path=outputs.stats_path
        end
        push!(validations, _validation_from_result_summary(data, outputs))
    end
    return validations
end

"""
    load_seed_msa(result; keepinserts=true, pid=nothing, index=nothing)

Load a selected ThorAxe seed MSA from an [`IdunaResult`](@ref) as a MIToS MSA.
"""
function load_seed_msa(result::Utils.IdunaResult; keepinserts::Bool = true,
        pid::Union{Nothing, Real} = nothing,
        index::Union{Nothing, Integer} = nothing)
    seed_index = _select_seed_index(result; pid, index, label = "seed")
    seed = result.thoraxe_msa.seeds[seed_index]
    seed_path = Utils._resolve_artifact_path(
        seed.stockholm_path, result.workdir)
    ResultsValidation.load_msa(seed_path; keepinserts)
end

"""
    load_expanded_msa(result; keepinserts=true, pid=nothing, index=nothing)

Load an expanded match-column MSA from an [`IdunaResult`](@ref) as a MIToS MSA.
"""
function load_expanded_msa(result::Utils.IdunaResult; keepinserts::Bool = true,
        pid::Union{Nothing, Real} = nothing,
        index::Union{Nothing, Integer} = nothing)
    isempty(result.expansions) &&
        error("This IdunaResult has no expanded MSA because it was run with no_expansion=true; use load_seed_msa(result) instead.")
    seed_index = _select_seed_index(result; pid, index, label = "expanded MSA")
    seed_index <= length(result.expansions) ||
        error("No expanded MSA is available at index $(seed_index).")
    expansion = result.expansions[seed_index]
    expansion === nothing &&
        error("No expanded MSA is available at index $(seed_index).")
    expansion_path = Utils._resolve_artifact_path(
        expansion.match_stockholm, result.workdir)
    ResultsValidation.load_msa(expansion_path; keepinserts)
end

function _select_seed_index(result::Utils.IdunaResult;
        pid::Union{Nothing, Real},
        index::Union{Nothing, Integer},
        label::AbstractString)
    (pid === nothing || index === nothing) ||
        error("Pass either pid or index to choose a $(label), not both.")
    seeds = result.thoraxe_msa.seeds
    isempty(seeds) && error("This IdunaResult has no selected seeds.")
    if index !== nothing
        1 <= index <= length(seeds) ||
            error("$(label) index $(index) is outside 1:$(length(seeds)).")
        return Int(index)
    end
    if pid !== nothing
        wanted = Float64(pid)
        found = findfirst(seed -> seed.pid == wanted, seeds)
        found === nothing && error("No selected $(label) is available for PID $(pid).")
        return found
    end
    length(seeds) == 1 ||
        error("This IdunaResult has $(length(seeds)) selected seeds; pass pid or index.")
    return 1
end
