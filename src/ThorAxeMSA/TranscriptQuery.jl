# Transcript Query Stage
# ----------------------

function _run_transcript_query_once(gene_core::AbstractString,
        workdir::AbstractString,
        species::Union{Nothing, AbstractString},
        specieslist::Union{Nothing, AbstractString},
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        orthology::AbstractString = "1:1",
        runner::Function = _thoraxe_runner(stdout_log, stderr_log))
    _orthology_relationships(orthology)
    cd(workdir) do
        if species === nothing
            ThorAxe.transcript_query(
                gene_core; orthology = orthology, specieslist = specieslist, runner = runner)
        else
            ThorAxe.transcript_query(
                gene_core; species = species, orthology = orthology,
                specieslist = specieslist, runner = runner)
        end
    end
    return nothing
end

function _ensure_cached_thoraxe_input(source_dir::AbstractString,
        workdir::AbstractString;
        overwrite::Bool = false,
        metadata)
    dest = _thoraxe_input_dir(workdir)
    source = abspath(String(source_dir))
    _has_valid_ensembl_bundle(source) ||
        error("Cached ThorAxe input at $(source) is not a valid transcript_query bundle.")
    source_fingerprint = _bundle_fingerprint(source)
    expected = merge(metadata,
        (;
            source_path = source,
            source_fingerprint
        ))
    identity = _thoraxe_input_identity(dest, expected)
    cache = _classify_thoraxe_input_stage(workdir, dest, identity; overwrite)
    reused = _maybe_reuse_thoraxe_input(workdir, dest, identity, expected, cache;
        overwrite,
        manifest_message = "Reusing manifest-backed ThorAxe transcript_query input.",
        legacy_message = "Adopting legacy ThorAxe transcript_query input.")
    reused === nothing || return reused

    _warn_stage_cache(cache, dest)
    action = _stage_action(cache)
    try
        _write_thoraxe_input_state(workdir, dest, :running, identity; action)
        if abspath(dest) != source
            isdir(dest) && safe_rm(dest, workdir)
            mkpath(dirname(dest))
            cp(source, dest; force = true)
        end
        _has_valid_ensembl_bundle(dest) ||
            error("Copied ThorAxe input bundle at $(dest) is incomplete.")
        _write_transcript_query_metadata!(dest, expected)
        _write_thoraxe_input_state(workdir, dest, :done,
            _thoraxe_input_identity(dest, expected); action)
    catch err
        _write_thoraxe_input_state(workdir, dest, :failed, identity;
            action,
            exception = _exception_summary(err))
        rethrow()
    end
    return dest
end

function _transcript_query_attempt_action!(gene_core::AbstractString,
        query_workdir::AbstractString,
        species::Union{Nothing, AbstractString},
        active_specieslist::Union{Nothing, AbstractString},
        stdout_log::AbstractString,
        stderr_log::AbstractString,
        tmp_gene_dir::AbstractString;
        attempt::Integer,
        attempts::Integer,
        orthology::AbstractString,
        runner::Function = _thoraxe_runner(stdout_log, stderr_log))
    start_details = (;
        gene = gene_core,
        species,
        orthology,
        attempt,
        attempts,
        _specieslist_log_details(active_specieslist)...)
    @debug "Starting ThorAxe transcript_query attempt." start_details...
    started = time()
    _run_transcript_query_once(
        gene_core, query_workdir, species, active_specieslist,
        stdout_log, stderr_log; orthology, runner)
    elapsed = round(time() - started; digits = 1)
    if _has_valid_ensembl_bundle(tmp_gene_dir)
        @debug "ThorAxe transcript_query attempt produced required input bundle." gene=gene_core attempt attempts elapsed_seconds=elapsed output_dir=tmp_gene_dir
        return :done
    end
    missing = _missing_transcript_query_outputs(tmp_gene_dir)
    invalid_details = (;
        gene = gene_core,
        attempt,
        attempts,
        elapsed_seconds = elapsed,
        missing_outputs = missing)
    if attempt >= attempts
        @warn "transcript_query produced an invalid bundle; no attempts remain." invalid_details...
        return :failed
    end
    @warn "transcript_query produced an invalid bundle; retrying with the same specieslist." invalid_details...
    return :retry
end

function _run_transcript_query_with_retries!(tmp_gene_dir::AbstractString,
        gene_core::AbstractString,
        query_workdir::AbstractString,
        species::Union{Nothing, AbstractString},
        initial_specieslist::Union{Nothing, AbstractString},
        attempts::Integer,
        stdout_log::AbstractString,
        stderr_log::AbstractString;
        orthology::AbstractString,
        runner::Function = _thoraxe_runner(stdout_log, stderr_log),
        sleep_fn::Function = sleep)
    for attempt in 1:attempts
        isdir(tmp_gene_dir) && rm(tmp_gene_dir; recursive = true, force = true)
        action = _transcript_query_attempt_action!(
            gene_core, query_workdir, species, initial_specieslist,
            stdout_log, stderr_log, tmp_gene_dir;
            attempt, attempts, orthology, runner)
        action === :done && return nothing
        action === :failed && break
        sleep_fn(_retry_wait_seconds(attempt))
    end
    return nothing
end

function _run_transcript_query_stage!(target::ResolvedTarget,
        workdir::AbstractString,
        input_dir::AbstractString,
        metadata;
        specieslist::Union{Nothing, AbstractString},
        max_retries::Integer,
        orthology::AbstractString,
        transcript_query_runner::Union{Nothing, Function},
        thoraxe_runner_factory::Function,
        sleep_fn::Function,
        action)
    isdir(input_dir) && safe_rm(input_dir, workdir)
    gene_core = strip_ensembl_version(target.ensembl_gene_id)
    logs = _thoraxe_logs_dir(workdir)
    stdout_log = joinpath(logs, "transcript_query_stdout.log")
    stderr_log = joinpath(logs, "transcript_query_stderr.log")
    species = _normalize_species_name(target.species)
    attempts = max(Int(max_retries), 1)
    active_specieslist = _normalized_specieslist(specieslist)
    run_details = (;
        species,
        attempts,
        logs_dir = logs)
    @info "Running ThorAxe transcript_query." run_details...
    stage_started = time()
    return mktempdir(workdir; prefix = "transcript_query_") do query_workdir
        tmp_gene_dir = joinpath(query_workdir, gene_core)
        runner = transcript_query_runner === nothing ?
                 thoraxe_runner_factory(stdout_log, stderr_log) : transcript_query_runner
        _run_transcript_query_with_retries!(
            tmp_gene_dir, gene_core, query_workdir, species, active_specieslist,
            attempts, stdout_log, stderr_log;
            orthology,
            runner,
            sleep_fn)

        _has_valid_ensembl_bundle(tmp_gene_dir) ||
            error("transcript_query did not create a valid Ensembl bundle at $(tmp_gene_dir). See $(stderr_log). If failures involve the species set, try a smaller curated specieslist.")
        mv(tmp_gene_dir, input_dir; force = true)
        _write_transcript_query_metadata!(input_dir, metadata)
        _write_thoraxe_input_state(workdir, input_dir, :done,
            _thoraxe_input_identity(input_dir, metadata); action)
        ready_details = (;
            elapsed_seconds = round(time() - stage_started; digits = 1))
        @info "ThorAxe transcript_query input is ready." ready_details...
        return input_dir
    end
end

function _ensure_transcript_query(target::ResolvedTarget, workdir::AbstractString;
        specieslist::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        cached_input_dir::Union{Nothing, AbstractString} = nothing,
        max_retries::Integer = 2,
        orthology::AbstractString = "1:1",
        transcript_query_runner::Union{Nothing, Function} = nothing,
        thoraxe_runner_factory::Function = _transcript_query_thoraxe_runner,
        sleep_fn::Function = sleep)
    _orthology_relationships(orthology)
    metadata = _expected_transcript_query_metadata(target;
        specieslist,
        orthology,
        source_kind = cached_input_dir === nothing ? "transcript_query" : "cached_input")
    if cached_input_dir !== nothing
        @info "Preparing ThorAxe input from cached transcript_query bundle." gene_id=target.ensembl_gene_id transcript_id=target.transcript_id source=cached_input_dir workdir
        return _ensure_cached_thoraxe_input(cached_input_dir, workdir; overwrite, metadata)
    end

    input_dir = _thoraxe_input_dir(workdir)
    identity = _thoraxe_input_identity(input_dir, metadata)
    cache = _classify_thoraxe_input_stage(workdir, input_dir, identity; overwrite)
    reused = _maybe_reuse_thoraxe_input(workdir, input_dir, identity, metadata, cache;
        overwrite,
        manifest_message = "Reusing manifest-backed ThorAxe transcript_query input.",
        legacy_message = "Reusing ThorAxe transcript_query input.")
    reused === nothing || return reused

    _warn_stage_cache(cache, input_dir)
    action = _stage_action(cache)
    try
        _write_thoraxe_input_state(workdir, input_dir, :running, identity; action)
        return _run_transcript_query_stage!(target, workdir, input_dir, metadata;
            specieslist,
            max_retries,
            orthology,
            transcript_query_runner,
            thoraxe_runner_factory,
            sleep_fn,
            action)
    catch err
        _write_thoraxe_input_state(workdir, input_dir, :failed, identity;
            action,
            exception = _exception_summary(err))
        rethrow()
    end
end

# Transcript Query Warnings
# -------------------------

function _species_from_biomart_errors_file(errors_path::AbstractString)
    isfile(errors_path) || return String[]
    species = String[]
    try
        for row in CSV.File(errors_path)
            value = getproperty(row, :Species)
            value === missing && continue
            normalized = _normalize_species_name(String(value))
            normalized === nothing && continue
            push!(species, normalized)
        end
    catch err
        err isa InterruptException && rethrow()
        return String[]
    end
    return _unique_nonempty_species(species)
end

function _species_from_biomart_stderr(stderr_log::AbstractString)
    isfile(stderr_log) || return String[]
    species = String[]
    for line in eachline(stderr_log)
        found = match(r"It can not found ([A-Za-z0-9_]+) in biomart", line)
        if found !== nothing
            push!(species, _normalize_species_name(found.captures[1]))
            continue
        end
        found = match(r"Download failed for \S+ in ([A-Za-z0-9_]+)", line)
        found !== nothing && push!(species, _normalize_species_name(found.captures[1]))
    end
    return _unique_nonempty_species(species)
end

function _biomart_transcript_query_warnings(input_dir::AbstractString,
        logs_dir::AbstractString)
    errors_path = joinpath(input_dir, "Ensembl", "errors.csv")
    stderr_log = joinpath(logs_dir, "transcript_query_stderr.log")
    error_species = _species_from_biomart_errors_file(errors_path)
    stderr_species = _species_from_biomart_stderr(stderr_log)
    species = _unique_nonempty_species(Iterators.flatten((error_species, stderr_species)))
    isempty(species) && return String[]
    return [
        "BioMart transcript_query failures recorded for species: $(join(species, ", ")). See $(errors_path) and $(stderr_log)."
    ]
end
