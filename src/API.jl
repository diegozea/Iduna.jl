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
datasets before `transcript_query`. That step can be slow, especially for broad
species lists; runs are often faster when a small curated `specieslist` is
provided. `pid_sample_count` controls how many PID-specific species samples are
drawn from each candidate `msa_0` and scored against that PID's full-species
MSA; the default is 45.
`pid_sample_fraction` controls the fraction of non-reference species retained in
each sample, and `pid_sample_seed` can make the sampling reproducible. If
omitted, a random seed is recorded with the result. Candidate `msa_0`
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
`no_expansion=true`.
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
        transcript_query_retries::Integer = 2,
        pid_sample_count::Integer = 45,
        pid_sample_fraction::Real = 0.8,
        pid_sample_seed::Union{Nothing, Integer} = nothing,
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
    validations = nothing
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
        Utils.write_json(joinpath(root, "target.json"), _target_summary(target, root))

        failed_stage = "thoraxe_msa"
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
            pid_sample_seed)

        expansions = Utils.ExpansionResult[]
        if !no_expansion
            failed_stage = "msa_expansion"
            expansions = [_expand_msa(target, seed, root;
                              mmseqs_db,
                              overwrite,
                              match_mode,
                              match_ratio,
                              hmmbuild_symfrac,
                              centroids,
                              threads)
                          for seed in thoraxe.seeds]
        end

        failed_stage = "validation"
        validations = [_validate_results(target, seed,
                           no_expansion ? nothing : expansions[index], root)
                       for (index, seed) in enumerate(thoraxe.seeds)]
        warnings = vcat(target.warnings, thoraxe.warnings,
            Iterators.flatten(validation.warnings for validation in validations)...)
        status = _pipeline_status(warnings)
        result = Utils.IdunaResult(;
            input_id = primary,
            workdir = root,
            target,
            thoraxe_msa = thoraxe,
            expansions,
            validations,
            warnings,
            status
        )
        result = Utils._relative_result_paths(result)
        Utils.write_json(joinpath(root, "result.json"), Utils.result_summary(result))
        return result
    catch err
        if root !== nothing
            _write_failure_result(
                primary, root, failed_stage, err; target, thoraxe, validations)
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
        validations = nothing)
    try
        Utils.write_json(joinpath(workdir, "result.json"),
            _failure_result_summary(input_id, workdir, failed_stage, err;
                target, thoraxe, validations))
    catch write_err
        @warn "Could not write failure result artifact." workdir exception=(write_err,
            catch_backtrace())
    end
    return nothing
end

function _failure_result_summary(input_id::AbstractString, workdir::AbstractString,
        failed_stage::AbstractString, err; target = nothing, thoraxe = nothing,
        validations = nothing)
    return (;
        input_id = String(input_id),
        workdir = String(workdir),
        status = "error",
        failed_stage = String(failed_stage),
        warnings = _partial_warnings(target, thoraxe, validations),
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
