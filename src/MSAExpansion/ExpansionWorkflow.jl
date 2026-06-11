# Expansion Workflow
# ------------------

function _expansion_context(target::ResolvedTarget,
        seed::SeedSelection,
        workdir::AbstractString,
        mmseqs_db::AbstractString;
        match_mode::Integer,
        match_ratio::Union{Nothing, Real},
        hmmbuild_symfrac::Real,
        centroids::Bool)
    seed_workdir = seed.workdir === nothing ? workdir : seed.workdir
    seed_stockholm = _resolve_artifact_path(seed.stockholm_path, seed_workdir)
    seed_fasta = seed.fasta_path === nothing ? nothing :
                 _resolve_artifact_path(seed.fasta_path, seed_workdir)
    run_dir = joinpath(_expansion_root(workdir), target.ensembl_gene_id,
        target.transcript_id, format_pid_dir(seed.pid))
    outputs = _expansion_output_paths(run_dir, target.transcript_id; centroids)
    identity = _expansion_identity(target, seed, seed_stockholm, seed_fasta, mmseqs_db;
        match_mode, match_ratio, hmmbuild_symfrac, centroids)
    return (;
        seed_stockholm,
        seed_fasta,
        run_dir,
        unpack_dir = joinpath(run_dir, "expanded_msa"),
        logs_dir = joinpath(run_dir, "logs"),
        outputs,
        identity,
        db_dir = joinpath(run_dir, "dbs"),
        seed_dir = joinpath(run_dir, "seeds"),
        hmm_dir = joinpath(run_dir, "hmm"),
        tmp_root = joinpath(run_dir, "tmp")
    )
end

function _prepare_expansion_cache!(run_dir::AbstractString,
        workdir::AbstractString,
        identity,
        outputs;
        overwrite::Bool)
    cache_warnings = String[]
    action = :run
    if overwrite && isdir(run_dir)
        @info "Clearing existing MSA expansion run directory." run_dir
        action = :rebuild
        safe_rm(run_dir, workdir)
    elseif !overwrite
        cache = _classify_step_state(run_dir, identity, outputs)
        if cache.reusable
            @info "Reusing cached MSA expansion." run_dir
            _write_step_state(run_dir, :done, identity, outputs; action = :reuse)
            return (; reusable = true, cache_warnings, action = :reuse)
        end
        action = cache.status === :missing ? :run : :rebuild
        if cache.warning !== nothing
            warning = String(cache.warning)
            push!(cache_warnings, warning)
            @warn warning run_dir status=cache.status
            isdir(run_dir) && safe_rm(run_dir, workdir)
        end
    end
    return (; reusable = false, cache_warnings, action)
end

function _cached_expansion_result(ctx, workdir::AbstractString)
    counts = _cached_hit_counts(ctx.outputs.hits_fasta, _seed_id_set(ctx.seed_stockholm))
    # The block table is derived from the saved alignments, so it can be restored lazily.
    _ensure_expansion_s_exon_blocks(ctx.outputs, ctx.identity.seed.pid)
    return ExpansionResult(;
        run_dir = ctx.run_dir,
        seed_stockholm = ctx.seed_stockholm,
        seed_fasta = ctx.seed_fasta,
        hits_fasta = ctx.outputs.hits_fasta,
        full_stockholm = ctx.outputs.full_stockholm,
        match_stockholm = ctx.outputs.match_stockholm,
        a3m_path = ctx.outputs.a3m_path,
        s_exon_blocks_tsv = ctx.outputs.s_exon_blocks_tsv,
        db_dir = ctx.db_dir,
        hmm_dir = ctx.hmm_dir,
        logs_dir = ctx.logs_dir,
        n_hits = counts.n_hits,
        n_new_hits = counts.n_new_hits,
        status = :skipped,
        workdir = String(workdir)
    )
end

function _prepare_expansion_dirs!(ctx, cache_warnings::Vector{String}, action)
    @debug "Preparing MSA expansion directories." run_dir=ctx.run_dir db_dir=ctx.db_dir hmm_dir=ctx.hmm_dir logs_dir=ctx.logs_dir
    mkpath.((
        ctx.db_dir, ctx.seed_dir, ctx.hmm_dir, ctx.tmp_root, ctx.unpack_dir, ctx.logs_dir))
    foreach(warning -> _write_cache_warning(ctx.logs_dir, warning), cache_warnings)
    _write_step_state(ctx.run_dir, :unfinished, ctx.identity, ctx.outputs;
        warnings = cache_warnings,
        action)
    return nothing
end

function _archive_expansion_seed(seed::SeedSelection, ctx)
    seed_label = "seed_pid$(format_pid(seed.pid))"
    run_dir = hasproperty(ctx, :run_dir) ? ctx.run_dir : nothing
    @debug "Archiving MSA expansion seed." pid=seed.pid seed_label run_dir
    archived_seed_sto = joinpath(ctx.seed_dir, "$(seed_label).sto")
    cp(ctx.seed_stockholm, archived_seed_sto; force = true)
    archived_seed_fasta = nothing
    if ctx.seed_fasta !== nothing && isfile(ctx.seed_fasta)
        archived_seed_fasta = joinpath(ctx.seed_dir, "$(seed_label).fasta")
        cp(ctx.seed_fasta, archived_seed_fasta; force = true)
    end
    # MMseqs cannot keep all Stockholm annotations, so archive them before conversion.
    sanitized_seed = prepare_stockholm_for_mmseqs(ctx.seed_stockholm,
        joinpath(ctx.seed_dir, "$(seed_label)_mmseqs.sto"))
    seed_alignment = read_file(archived_seed_sto, Stockholm; keepinserts = true)
    seed_names = String.(sequencenames(seed_alignment))
    # Save the seed's s-exon labels before temporary MMseqs files strip them away.
    seed_s_exon_codes = s_exon_codes(seed_alignment)
    seed_s_exon_code_map = sort!(collect(s_exon_code_map(seed_alignment)); by = first)
    return (;
        archived_seed_sto,
        archived_seed_fasta,
        sanitized_seed,
        seed_names,
        seed_pid = seed.pid,
        seed_s_exon_codes,
        seed_s_exon_code_map,
        seed_set = Set(_normalize_id.(seed_names))
    )
end

function _write_expansion_hits!(hits_tsv::AbstractString,
        run_dir::AbstractString,
        hmm_dir::AbstractString,
        seed_set::Set{String})
    all_hits, filtered_hits = _collect_hits(hits_tsv, seed_set)
    raw_hits_fasta = joinpath(run_dir, "mmseqs_hits_raw.fasta")
    filtered_fasta = joinpath(hmm_dir, "mmseqs_hits_filtered.fasta")
    write_fasta(raw_hits_fasta, all_hits)
    write_fasta(filtered_fasta, filtered_hits)
    return (; all_hits, filtered_hits, raw_hits_fasta, filtered_fasta)
end

function _build_seed_hmm!(seed_dir::AbstractString,
        sanitized_seed::AbstractString,
        hmmbuild_symfrac::Real,
        logs_dir::AbstractString)
    hmm_path = joinpath(seed_dir, "seed.hmm")
    annotated_seed = joinpath(seed_dir, "seed_annotated.sto")
    _run_labeled(
        `$(HMMER_jll.hmmbuild()) --symfrac $(Float64(hmmbuild_symfrac)) -O $annotated_seed $hmm_path $sanitized_seed`,
        "hmmbuild", logs_dir)
    normalize_stockholm_annotations!(annotated_seed)
    return (; hmm_path, annotated_seed)
end

function _align_expansion_hits!(hmm_dir::AbstractString,
        annotated_seed::AbstractString,
        sanitized_seed::AbstractString,
        hmm_path::AbstractString,
        filtered_fasta::AbstractString,
        filtered_hits,
        seed_set::Set{String},
        logs_dir::AbstractString)
    aligned_sto = joinpath(hmm_dir, "alignment_with_hits.sto")
    hits_aligned_sto = joinpath(hmm_dir, "aligned_hits_only.sto")
    if isempty(filtered_hits)
        # With no new hits, the expanded alignment is just the annotated seed.
        cp(annotated_seed, aligned_sto; force = true)
    else
        _run_labeled(
            `$(HMMER_jll.hmmalign()) --mapali $sanitized_seed --trim --outformat stockholm -o $aligned_sto $hmm_path $filtered_fasta`,
            "hmmalign", logs_dir)
        aligned_msa = read_file(aligned_sto, Stockholm; keepinserts = true)
        hit_indices = [i
                       for (i, name) in enumerate(sequencenames(aligned_msa))
                       if !(_normalize_id(String(name)) in seed_set)]
        write_file(hits_aligned_sto, aligned_msa[hit_indices, :], Stockholm)
    end
    normalize_stockholm_annotations!(aligned_sto)
    return aligned_sto
end

function _write_expansion_alignment_outputs!(target::ResolvedTarget,
        unpack_dir::AbstractString,
        aligned_sto::AbstractString,
        seed_names::Vector{String},
        raw_hits_fasta::AbstractString,
        archived,
        pid::Real)
    match_stockholm = joinpath(unpack_dir, "$(target.transcript_id)_matchonly.sto")
    open(match_stockholm, "w") do io
        run(pipeline(`$(HMMER_jll.esl_alimask()) --rf-is-mask $aligned_sto`, stdout = io))
    end
    normalize_stockholm_annotations!(match_stockholm)

    full_stockholm = joinpath(unpack_dir, "$(target.transcript_id)_full.sto")
    full_alignment = _reorder_alignment(
        read_file(aligned_sto, Stockholm; keepinserts = true), seed_names)
    match_alignment = _reorder_alignment(
        read_file(match_stockholm, Stockholm; keepinserts = true), seed_names)
    _restore_s_exon_annotations!(full_alignment, archived)
    _restore_s_exon_annotations!(match_alignment, archived)
    write_file(full_stockholm, full_alignment, Stockholm)
    write_file(match_stockholm, match_alignment, Stockholm)

    a3m_path = joinpath(unpack_dir, "$(target.transcript_id)_expanded.a3m")
    write_file(a3m_path, match_alignment, A3M)
    s_exon_blocks_tsv = joinpath(unpack_dir, "$(target.transcript_id)_s_exon_blocks.tsv")
    _write_expansion_s_exon_blocks(
        s_exon_blocks_tsv, full_alignment, match_alignment, pid)
    hits_copy = joinpath(unpack_dir, "$(target.transcript_id)_hits_raw.fasta")
    cp(raw_hits_fasta, hits_copy; force = true)
    return (; match_stockholm, full_stockholm, a3m_path, s_exon_blocks_tsv, hits_copy)
end

function _write_centroids_if_requested!(target::ResolvedTarget,
        mmseqs_db::AbstractString,
        db_paths,
        tmp_dir::AbstractString,
        ctx,
        archived,
        hmm_paths,
        centroids::Bool)
    centroids || return nothing
    _write_centroid_msa(target.transcript_id,
        mmseqs_db,
        db_paths,
        tmp_dir,
        joinpath(ctx.run_dir, "centroid_msa"),
        archived.seed_set,
        archived.seed_names,
        archived.sanitized_seed,
        hmm_paths.hmm_path,
        hmm_paths.annotated_seed,
        archived,
        ctx.logs_dir)
    return nothing
end

function _run_expansion_workflow!(target::ResolvedTarget,
        ctx,
        archived,
        mmseqs_db::AbstractString;
        match_mode::Integer,
        match_ratio::Union{Nothing, Real},
        hmmbuild_symfrac::Real,
        centroids::Bool,
        mmseqs_threads::Union{Nothing, Integer})
    return mktempdir(ctx.tmp_root; prefix = "mmseqs_tmp_") do tmp_dir
        @info "Running MMseqs expansion workflow." mmseqs_threads
        @debug "Using MMseqs expansion temporary directory." pid=archived.seed_pid tmp_dir
        # MMseqs finds candidate homologs; HMMER maps them back to seed columns.
        db_paths = _mmseqs_search(archived.sanitized_seed, mmseqs_db, ctx.db_dir,
            tmp_dir, ctx.logs_dir; match_mode, match_ratio, mmseqs_threads)
        hits_tsv = joinpath(ctx.run_dir, "mmseqs_hits_raw.tsv")
        @debug "Converting MMseqs hits." pid=archived.seed_pid hits_tsv
        _run_labeled(
            `$(MMseqs2_jll.mmseqs()) convertalis $(db_paths.profile_db) $(db_paths.seq_db) $(db_paths.realigned_result_db) $hits_tsv --format-output query,target,tseq`,
            "convertalis", ctx.logs_dir)

        @debug "Writing MSA expansion hits." pid=archived.seed_pid
        hits = _write_expansion_hits!(
            hits_tsv, ctx.run_dir, ctx.hmm_dir, archived.seed_set)
        @debug "Building seed HMM." pid=archived.seed_pid hmmbuild_symfrac
        hmm_paths = _build_seed_hmm!(
            ctx.seed_dir, archived.sanitized_seed, hmmbuild_symfrac, ctx.logs_dir)
        archived = _with_seed_match_s_exon_codes(archived, hmm_paths.annotated_seed)
        @debug "Aligning MSA expansion hits." pid=archived.seed_pid n_new_hits=length(hits.filtered_hits)
        aligned_sto = _align_expansion_hits!(ctx.hmm_dir,
            hmm_paths.annotated_seed,
            archived.sanitized_seed,
            hmm_paths.hmm_path,
            hits.filtered_fasta,
            hits.filtered_hits,
            archived.seed_set,
            ctx.logs_dir)
        @debug "Writing expanded alignment outputs." pid=archived.seed_pid
        output_paths = _write_expansion_alignment_outputs!(
            target, ctx.unpack_dir, aligned_sto, archived.seed_names,
            hits.raw_hits_fasta, archived, archived.seed_pid)
        _write_centroids_if_requested!(
            target, mmseqs_db, db_paths, tmp_dir, ctx, archived, hmm_paths, centroids)
        @info "MSA expansion completed." pid=archived.seed_pid n_hits=length(hits.all_hits) n_new_hits=length(hits.filtered_hits)
        return (;
            output_paths,
            n_hits = length(hits.all_hits),
            n_new_hits = length(hits.filtered_hits)
        )
    end
end

function _finished_expansion_result(ctx,
        archived,
        run_outputs,
        workdir::AbstractString)
    return ExpansionResult(;
        run_dir = ctx.run_dir,
        seed_stockholm = archived.archived_seed_sto,
        seed_fasta = archived.archived_seed_fasta,
        hits_fasta = run_outputs.output_paths.hits_copy,
        full_stockholm = run_outputs.output_paths.full_stockholm,
        match_stockholm = run_outputs.output_paths.match_stockholm,
        a3m_path = run_outputs.output_paths.a3m_path,
        s_exon_blocks_tsv = run_outputs.output_paths.s_exon_blocks_tsv,
        db_dir = ctx.db_dir,
        hmm_dir = ctx.hmm_dir,
        logs_dir = ctx.logs_dir,
        n_hits = run_outputs.n_hits,
        n_new_hits = run_outputs.n_new_hits,
        status = :ok,
        workdir = String(workdir)
    )
end
