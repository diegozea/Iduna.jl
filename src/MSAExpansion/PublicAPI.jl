# Public API
# ----------

"""
    expand_msa(target, seed, workdir; mmseqs_db, mmseqs_db_manifest=nothing,
               overwrite=false, match_mode=1,
               match_ratio=nothing, hmmbuild_symfrac=0.0, centroids=false,
               mmseqs_threads=Threads.nthreads()) -> ExpansionResult

Expand one selected seed MSA against an MMseqs2 database.

# Arguments

- `target::ResolvedTarget`: resolved target metadata.
- `seed::SeedSelection`: selected ThorAxe seed to expand.
- `workdir::AbstractString`: Iduna work directory.

# Keywords

- `mmseqs_db::AbstractString`: MMseqs2 database prefix.
- `mmseqs_db_manifest = nothing`: manifest declaring an immutable database
  snapshot and SHA-256 checksums of all its component files. Without a valid
  manifest, expansion runs but cached expansion results are not reused. Iduna
  trusts the snapshot's immutability instead of rehashing the database each run.
- `overwrite::Bool = false`: reuse package-owned expansion outputs when their
  implementation, tools, database snapshot, thread settings, run identity and
  required-file digests still match. When `true`, Iduna rebuilds those outputs.
  Older states without this provenance are rebuilt on the next run.
- `match_mode::Integer = 1`: MMseqs2 profile matching mode passed to
  `mmseqs msa2profile --match-mode`.
- `match_ratio = nothing`: optional MMseqs2 match-ratio setting. When `nothing`,
  Iduna does not pass `--match-ratio` to MMseqs2. Other values are passed to
  `mmseqs msa2profile --match-ratio`.
- `hmmbuild_symfrac::Real = 0.0`: HMMER `hmmbuild` symbol fraction. Must be
  between `0.0` and `1.0`.
- `centroids::Bool = false`: only save the full expanded MSA. When `true`, Iduna
  also saves a centroid-level MSA.
- `mmseqs_threads = Threads.nthreads()`: number of threads passed to MMseqs2.
  This uses the active Julia thread count unless another value is passed.
  Passing `nothing` leaves the tool's thread default unspecified and disables reuse.

# Returns

- An [`ExpansionResult`](@ref) with paths to the expanded MSA and logs.
"""
function expand_msa(target::ResolvedTarget,
        seed::SeedSelection,
        workdir::AbstractString;
        mmseqs_db::AbstractString,
        mmseqs_db_manifest::Union{Nothing, AbstractString} = nothing,
        overwrite::Bool = false,
        match_mode::Integer = 1,
        match_ratio::Union{Nothing, Real} = nothing,
        hmmbuild_symfrac::Real = 0.0,
        centroids::Bool = false,
        mmseqs_threads::Union{Nothing, Integer} = Threads.nthreads())
    0.0 <= hmmbuild_symfrac <= 1.0 ||
        error("hmmbuild_symfrac must be between 0.0 and 1.0.")
    @info "Preparing MSA expansion." pid=seed.pid mmseqs_db centroids
    ensure_mmseqs_db(mmseqs_db)

    ctx = _expansion_context(target, seed, workdir, mmseqs_db;
        match_mode, match_ratio, hmmbuild_symfrac, centroids,
        mmseqs_db_manifest, mmseqs_threads)
    cache = _prepare_expansion_cache!(ctx.run_dir, workdir, ctx.identity, ctx.outputs;
        overwrite)
    cache.reusable && return _cached_expansion_result(ctx, workdir)

    _prepare_expansion_dirs!(ctx, cache.cache_warnings, cache.action)
    archived = _archive_expansion_seed(seed, ctx)

    try
        run_outputs = _run_expansion_workflow!(target, ctx, archived, mmseqs_db;
            match_mode, match_ratio, hmmbuild_symfrac, centroids, mmseqs_threads)
        _write_step_state(ctx.run_dir, :done, ctx.identity, ctx.outputs;
            warnings = cache.cache_warnings,
            action = cache.action)
        @info "MSA expansion completed." pid=archived.seed_pid n_hits=run_outputs.n_hits n_new_hits=run_outputs.n_new_hits
        return _finished_expansion_result(ctx, archived, run_outputs, workdir)
    catch err
        err isa InterruptException && rethrow()
        _write_step_state(ctx.run_dir, :failed, ctx.identity, ctx.outputs;
            warnings = cache.cache_warnings,
            exception = _exception_summary(err),
            action = cache.action)
        rethrow()
    end
end
