using Base.Threads

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
`"m:n"`). By default, Iduna filters the ThorAxe species list with Ensembl
homology and then filters against currently available BioMart Ensembl Gene
datasets before `transcript_query`. `transcript_query_timeout_seconds` limits
each Ensembl query attempt. Iduna retries that step and can retry without
`specieslist` after a timeout. Set `thoraxe_timeout_seconds` to give the
baseline and PID ThorAxe runs the same kind of wall-clock guard. Pass
`thoraxe_input_dir` to reuse a complete transcript_query bundle instead of
fetching it again. Set `centroids=true` to also save the centroid-level MSA
before MMseqs2 expands centroid hits to cluster members; the main expansion and
validation still use the full expanded MSA. `centroids=true` requires
expansion and cannot be combined with `no_expansion=true`.
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
        specieslist::Union{Nothing, AbstractString} = nothing,
        orthology::AbstractString = "1:1",
        specieslist_filter::Bool = true,
        biomart_datasets_filter::Bool = true,
        thoraxe_input_dir::Union{Nothing, AbstractString} = nothing,
        transcript_query_timeout_seconds::Union{Nothing, Real} = 180,
        transcript_query_timeout_max_seconds::Union{Nothing, Real} = 240,
        transcript_query_retries::Integer = 2,
        allow_specieslist_timeout_fallback::Bool = true,
        thoraxe_timeout_seconds::Union{Nothing, Real} = nothing,
        match_mode::Integer = 1,
        match_ratio::Union{Nothing, Real} = nothing,
        hmmbuild_symfrac::Real = 0.0,
        centroids::Bool = false,
        threads::Union{Nothing, Integer} = Threads.nthreads(),
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
    validation = nothing
    failed_stage = "prepare_output_dir"

    try
        root = Utils.prepare_output_dir(primary; workdir, output_dir, overwrite)

        failed_stage = "resolve_target"
        target = _resolve_target(primary;
            workdir = root,
            uniprot_id = supplied_uniprot,
            ensembl_gene_id,
            ensembl_protein_id,
            transcript_id = disambiguating_transcript,
            species)
        Utils.write_json(joinpath(root, "target.json"), _target_summary(target))

        failed_stage = "thoraxe_msa"
        thoraxe = _build_thoraxe_msa(target, root;
            pid_thresholds,
            specieslist,
            orthology,
            specieslist_filter,
            biomart_datasets_filter,
            cached_thoraxe_input_dir = thoraxe_input_dir,
            overwrite,
            transcript_query_timeout_seconds,
            transcript_query_timeout_max_seconds,
            transcript_query_retries,
            allow_specieslist_timeout_fallback,
            thoraxe_timeout_seconds)

        expansion = nothing
        if !no_expansion
            failed_stage = "msa_expansion"
            expansion = _expand_msa(target, thoraxe.best_seed, root;
                mmseqs_db,
                overwrite,
                match_mode,
                match_ratio,
                hmmbuild_symfrac,
                centroids,
                threads)
        end

        failed_stage = "validation"
        validation = _validate_results(target, thoraxe.best_seed, expansion, root)
        warnings = vcat(target.warnings, thoraxe.warnings, validation.warnings)
        status = _pipeline_status(warnings)
        result = Utils.IdunaResult(;
            input_id = primary,
            workdir = root,
            target,
            thoraxe_msa = thoraxe,
            expansion,
            validation,
            warnings,
            status
        )
        Utils.write_json(joinpath(root, "result.json"), Utils.result_summary(result))
        return result
    catch err
        if root !== nothing
            _write_failure_result(
                primary, root, failed_stage, err; target, thoraxe, validation)
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

function _write_failure_result(input_id::AbstractString, workdir::AbstractString,
        failed_stage::AbstractString, err; target = nothing, thoraxe = nothing,
        validation = nothing)
    try
        Utils.write_json(joinpath(workdir, "result.json"),
            _failure_result_summary(input_id, workdir, failed_stage, err;
                target, thoraxe, validation))
    catch write_err
        @warn "Could not write failure result artifact." workdir exception=(write_err,
            catch_backtrace())
    end
    return nothing
end

function _failure_result_summary(input_id::AbstractString, workdir::AbstractString,
        failed_stage::AbstractString, err; target = nothing, thoraxe = nothing,
        validation = nothing)
    return (;
        input_id = String(input_id),
        workdir = String(workdir),
        status = "error",
        failed_stage = String(failed_stage),
        warnings = _partial_warnings(target, thoraxe, validation),
        target = target === nothing ? nothing : _target_summary(target),
        exception = _exception_summary(err)
    )
end

function _partial_warnings(target, thoraxe, validation)
    warnings = String[]
    target !== nothing && append!(warnings, target.warnings)
    thoraxe !== nothing && append!(warnings, thoraxe.warnings)
    validation !== nothing && append!(warnings, validation.warnings)
    return warnings
end

function _exception_summary(err)
    summary = (;
        type = string(typeof(err)),
        message = sprint(showerror, err)
    )
    if err isa ThorAxeMSA._CommandTimeoutError
        return merge(
            summary, (;
                command = err.command,
                stdout_log = err.stdout_log,
                stderr_log = err.stderr_log
            ))
    end
    return summary
end

function _normalize_primary_input(; id, uniprot_id, ensembl_transcript_id, transcript_id)
    # Accept old and new argument names, but still require one clear input.
    provided_primary = Pair{Symbol, String}[]
    id !== nothing && push!(provided_primary, :id => String(id))
    ensembl_transcript_id !== nothing &&
        push!(provided_primary, :ensembl_transcript_id => String(ensembl_transcript_id))
    if isempty(provided_primary) && uniprot_id !== nothing
        push!(provided_primary, :uniprot_id => String(uniprot_id))
    end

    if isempty(provided_primary)
        transcript_id === nothing &&
            error("Pass id, uniprot_id, ensembl_transcript_id, or transcript_id.")
        return String(transcript_id), nothing, nothing
    end
    first_value = first(provided_primary).second
    if any(pair -> pair.second != first_value, provided_primary)
        error("Conflicting primary identifiers were provided: $(provided_primary).")
    end

    kind = Utils.id_kind(first_value)
    # For UniProt input, transcript_id is only a tie-breaker among candidates.
    if kind === :uniprot && uniprot_id !== nothing && String(uniprot_id) != first_value
        error("Conflicting UniProt identifiers were provided: id=$(first_value), uniprot_id=$(uniprot_id).")
    end
    disambiguating_transcript = kind === :uniprot ? transcript_id : nothing
    supplied_uniprot = kind === :uniprot ? first_value :
                       (uniprot_id === nothing ? nothing : String(uniprot_id))
    return first_value, disambiguating_transcript, supplied_uniprot
end

function _target_summary(target::Utils.ResolvedTarget)
    return (;
        input_id = target.input_id,
        input_kind = String(target.input_kind),
        uniprot_id = target.uniprot_id,
        ensembl_gene_id = target.ensembl_gene_id,
        transcript_id = target.transcript_id,
        ensembl_protein_id = target.ensembl_protein_id,
        species = target.species,
        uniprot_sequence_path = target.uniprot_sequence_path,
        ensembl_protein_sequence_path = target.ensembl_protein_sequence_path,
        sequence_validated = target.sequence_validated,
        mapping_confirmed = target.mapping_confirmed,
        warnings = target.warnings
    )
end

"""
    load_seed_msa(result; keepinserts=true)

Load the selected ThorAxe seed MSA from an [`IdunaResult`](@ref) as a MIToS MSA.
"""
function load_seed_msa(result::Utils.IdunaResult; keepinserts::Bool = true)
    ResultsValidation.load_seed_msa(result.thoraxe_msa.best_seed; keepinserts)
end

"""
    load_expanded_msa(result; keepinserts=true)

Load the expanded match-column MSA from an [`IdunaResult`](@ref) as a MIToS MSA.
"""
function load_expanded_msa(result::Utils.IdunaResult; keepinserts::Bool = true)
    result.expansion === nothing &&
        error("This IdunaResult has no expanded MSA because it was run with no_expansion=true; use load_seed_msa(result) instead.")
    ResultsValidation.load_expanded_msa(result.expansion; keepinserts)
end
