module MSAExpansion

using HMMER_jll
using MIToS.MSA
using MMseqs2_jll

using ..Utils: ExpansionResult, ResolvedTarget, SeedSelection, ensure_mmseqs_db,
               format_pid, run_logged, safe_rm, write_fasta

export expand_msa,
       normalize_stockholm_annotations!,
       prepare_stockholm_for_mmseqs

_normalize_id(s::AbstractString) = String(split(String(s))[1])

function _expansion_root(workdir::AbstractString)
    joinpath(workdir, "expansion")
end

function _expanded_outputs_exist(dir::AbstractString, transcript_id::AbstractString)
    all(isfile,
        (
            joinpath(dir, "$(transcript_id)_full.sto"),
            joinpath(dir, "$(transcript_id)_matchonly.sto"),
            joinpath(dir, "$(transcript_id)_expanded.a3m"),
            joinpath(dir, "$(transcript_id)_hits_raw.fasta")
        ))
end

function prepare_stockholm_for_mmseqs(source::AbstractString, dest::AbstractString)
    mkpath(dirname(dest))
    lines = readlines(source)
    has_header = !isempty(lines) && startswith(strip(lines[1]), "# STOCKHOLM")
    open(dest, "w") do io
        has_header || println(io, "# STOCKHOLM 1.0")
        for line in lines
            println(io, line)
        end
        if isempty(lines) || strip(lines[end]) != "//"
            println(io, "//")
        end
    end
    return dest
end

function normalize_stockholm_annotations!(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return path

    comments = String[]
    gf_order = String[]
    gf_data = Dict{String, Vector{String}}()
    gs_order = Tuple{String, String}[]
    gs_data = Dict{Tuple{String, String}, Vector{String}}()
    gc_order = String[]
    gc_data = Dict{String, String}()
    gr_order = Tuple{String, String}[]
    gr_data = Dict{Tuple{String, String}, String}()
    seq_order = String[]
    seq_data = Dict{String, String}()

    for line in lines
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, "# STOCKHOLM") || startswith(stripped, "//")
            continue
        elseif startswith(line, "#=GF")
            parts = split(line; limit = 3)
            length(parts) < 3 && continue
            feature = parts[2]
            haskey(gf_data, feature) ||
                (push!(gf_order, feature); gf_data[feature] = String[])
            push!(gf_data[feature], parts[3])
        elseif startswith(line, "#=GS")
            parts = split(line; limit = 4)
            length(parts) < 4 && continue
            key = (parts[2], parts[3])
            haskey(gs_data, key) || (push!(gs_order, key); gs_data[key] = String[])
            push!(gs_data[key], parts[4])
        elseif startswith(line, "#=GC")
            parts = split(line; limit = 3)
            length(parts) < 3 && continue
            feature = parts[2]
            data = replace(parts[3], ' ' => "")
            if haskey(gc_data, feature)
                gc_data[feature] = string(gc_data[feature], data)
            else
                push!(gc_order, feature)
                gc_data[feature] = data
            end
        elseif startswith(line, "#=GR")
            parts = split(line; limit = 4)
            length(parts) < 4 && continue
            key = (parts[2], parts[3])
            data = replace(parts[4], ' ' => "")
            if haskey(gr_data, key)
                gr_data[key] = string(gr_data[key], data)
            else
                push!(gr_order, key)
                gr_data[key] = data
            end
        elseif startswith(line, '#')
            push!(comments, stripped)
        else
            parts = split(line; limit = 2)
            length(parts) < 2 && continue
            name = strip(parts[1])
            fragment = replace(strip(parts[2]), ' ' => "")
            if haskey(seq_data, name)
                seq_data[name] = string(seq_data[name], fragment)
            else
                push!(seq_order, name)
                seq_data[name] = fragment
            end
        end
    end

    open(path, "w") do io
        println(io, "# STOCKHOLM 1.0")
        foreach(line -> println(io, line), comments)
        for feature in gf_order, value in gf_data[feature]

            println(io, "#=GF ", feature, ' ', value)
        end
        for key in gs_order, value in gs_data[key]

            println(io, "#=GS ", key[1], ' ', key[2], ' ', value)
        end
        for feature in gc_order
            println(io, "#=GC ", feature, ' ', gc_data[feature])
        end
        for name in seq_order
            println(io, name, '\t', seq_data[name])
            for key in gr_order
                key[1] == name || continue
                println(io, "#=GR ", key[1], ' ', key[2], ' ', gr_data[key])
            end
        end
        println(io, "//")
    end
    return path
end

function _run_labeled(cmd::Cmd, label::AbstractString, logs_dir::AbstractString)
    mkpath(logs_dir)
    run_logged(cmd;
        stdout_path = joinpath(logs_dir, "$(label)_stdout.log"),
        stderr_path = joinpath(logs_dir, "$(label)_stderr.log"))
    return nothing
end

function _mmseqs_search(seed_msa_sto::AbstractString,
        base_db::AbstractString,
        db_dir::AbstractString,
        tmp_dir::AbstractString,
        logs_dir::AbstractString;
        match_mode::Integer = 1,
        match_ratio::Union{Nothing, Real} = nothing,
        threads::Union{Nothing, Integer} = nothing)
    mmseqs_bin = MMseqs2_jll.mmseqs()
    mkpath(db_dir)
    seed_db = joinpath(db_dir, "seed_msa_db")
    profile_db = joinpath(db_dir, "seed_profile_db")
    search_result_db = joinpath(db_dir, "search_results_db")
    expanded_result_db = joinpath(db_dir, "expanded_results_db")
    realigned_result_db = joinpath(db_dir, "realigned_results_db")
    aln_db = string(base_db, "_aln")
    seq_db = string(base_db, "_seq")

    _run_labeled(`$(mmseqs_bin) convertmsa $seed_msa_sto $seed_db`, "convertmsa", logs_dir)

    profile_cmd = `$(mmseqs_bin) msa2profile $seed_db $profile_db --match-mode $(match_mode)`
    match_ratio !== nothing &&
        (profile_cmd = `$profile_cmd --match-ratio $(Float64(match_ratio))`)
    _run_labeled(profile_cmd, "msa2profile", logs_dir)

    search_cmd = `$(mmseqs_bin) search $profile_db $base_db $search_result_db $tmp_dir -a`
    threads !== nothing && (search_cmd = `$search_cmd --threads $(Int(threads))`)
    _run_labeled(search_cmd, "search", logs_dir)

    expandaln_cmd = `$(mmseqs_bin) expandaln $profile_db $base_db $search_result_db $aln_db $expanded_result_db`
    threads !== nothing && (expandaln_cmd = `$expandaln_cmd --threads $(Int(threads))`)
    _run_labeled(expandaln_cmd, "expandaln", logs_dir)

    align_cmd = `$(mmseqs_bin) align $profile_db $seq_db $expanded_result_db $realigned_result_db`
    threads !== nothing && (align_cmd = `$align_cmd --threads $(Int(threads))`)
    _run_labeled(align_cmd, "align", logs_dir)

    return (; seed_db, profile_db, seq_db, search_result_db,
        expanded_result_db, realigned_result_db)
end

function _reorder_alignment(msa::AbstractMultipleSequenceAlignment, priority::Vector{String})
    name_to_idx = Dict(String(name) => idx for (idx, name) in enumerate(sequencenames(msa)))
    priority_idx = Int[]
    seen = Set{Int}()
    for name in priority
        idx = get(name_to_idx, String(name), nothing)
        idx === nothing && continue
        push!(priority_idx, idx)
        push!(seen, idx)
    end
    remaining = [idx for idx in 1:nsequences(msa) if !(idx in seen)]
    return msa[vcat(priority_idx, remaining), :]
end

function _collect_hits(hits_tsv::AbstractString, seed_set::Set{String})
    all_hits = Tuple{String, String}[]
    filtered_hits = Tuple{String, String}[]
    seen_all = Set{String}()
    seen_filtered = Set{String}()
    open(hits_tsv, "r") do io
        for line in eachline(io)
            parts = split(line, '\t')
            length(parts) < 3 && continue
            target_name = _normalize_id(strip(parts[2]))
            seq = uppercase(replace(strip(parts[3]), '-' => ""))
            isempty(target_name) && continue
            occursin(r"[A-Za-z]", seq) || continue
            if !(target_name in seen_all)
                push!(all_hits, (target_name, seq))
                push!(seen_all, target_name)
            end
            if target_name in seed_set || target_name in seen_filtered
                continue
            end
            push!(filtered_hits, (target_name, seq))
            push!(seen_filtered, target_name)
        end
    end
    return all_hits, filtered_hits
end

function _seed_id_set(seed_stockholm::AbstractString)
    seed_alignment = read_file(seed_stockholm, Stockholm; keepinserts = true)
    return Set(_normalize_id.(String.(sequencenames(seed_alignment))))
end

function _cached_hit_counts(hits_fasta::AbstractString, seed_set::Set{String})
    if filesize(hits_fasta) == 0
        return (; n_hits = 0, n_new_hits = 0)
    end
    hits_msa = read_file(hits_fasta, FASTA)
    hit_names = _normalize_id.(String.(sequencenames(hits_msa)))
    return (;
        n_hits = nsequences(hits_msa),
        n_new_hits = count(name -> !(name in seed_set), hit_names)
    )
end

function expand_msa(target::ResolvedTarget,
        seed::SeedSelection,
        workdir::AbstractString;
        mmseqs_db::AbstractString,
        overwrite::Bool = false,
        match_mode::Integer = 1,
        match_ratio::Union{Nothing, Real} = nothing,
        hmmbuild_symfrac::Real = 0.0,
        threads::Union{Nothing, Integer} = Threads.nthreads())
    ensure_mmseqs_db(mmseqs_db)

    run_dir = joinpath(_expansion_root(workdir), target.ensembl_gene_id, target.transcript_id)
    unpack_dir = joinpath(run_dir, "expanded_msa")
    logs_dir = joinpath(run_dir, "logs")
    if overwrite && isdir(run_dir)
        safe_rm(run_dir, workdir)
    elseif !overwrite && _expanded_outputs_exist(unpack_dir, target.transcript_id)
        hits_fasta = joinpath(unpack_dir, "$(target.transcript_id)_hits_raw.fasta")
        counts = _cached_hit_counts(hits_fasta, _seed_id_set(seed.stockholm_path))
        return ExpansionResult(;
            run_dir,
            seed_stockholm = seed.stockholm_path,
            seed_fasta = seed.fasta_path,
            hits_fasta,
            full_stockholm = joinpath(unpack_dir, "$(target.transcript_id)_full.sto"),
            match_stockholm = joinpath(unpack_dir, "$(target.transcript_id)_matchonly.sto"),
            a3m_path = joinpath(unpack_dir, "$(target.transcript_id)_expanded.a3m"),
            db_dir = joinpath(run_dir, "dbs"),
            hmm_dir = joinpath(run_dir, "hmm"),
            logs_dir,
            n_hits = counts.n_hits,
            n_new_hits = counts.n_new_hits,
            status = :skipped
        )
    end

    db_dir = joinpath(run_dir, "dbs")
    seed_dir = joinpath(run_dir, "seeds")
    hmm_dir = joinpath(run_dir, "hmm")
    tmp_root = joinpath(run_dir, "tmp")
    mkpath.((db_dir, seed_dir, hmm_dir, tmp_root, unpack_dir, logs_dir))

    seed_label = "seed_pid$(format_pid(seed.pid))"
    archived_seed_sto = joinpath(seed_dir, "$(seed_label).sto")
    cp(seed.stockholm_path, archived_seed_sto; force = true)
    archived_seed_fasta = nothing
    if seed.fasta_path !== nothing && isfile(seed.fasta_path)
        archived_seed_fasta = joinpath(seed_dir, "$(seed_label).fasta")
        cp(seed.fasta_path, archived_seed_fasta; force = true)
    end
    sanitized_seed = prepare_stockholm_for_mmseqs(seed.stockholm_path,
        joinpath(seed_dir, "$(seed_label)_mmseqs.sto"))
    seed_alignment = read_file(archived_seed_sto, Stockholm; keepinserts = true)
    seed_names = String.(sequencenames(seed_alignment))
    seed_set = Set(_normalize_id.(seed_names))

    tmp_dir = mktempdir(tmp_root; prefix = "mmseqs_tmp_")
    try
        db_paths = _mmseqs_search(sanitized_seed, mmseqs_db, db_dir, tmp_dir, logs_dir;
            match_mode, match_ratio, threads)
        hits_tsv = joinpath(run_dir, "mmseqs_hits_raw.tsv")
        _run_labeled(
            `$(MMseqs2_jll.mmseqs()) convertalis $(db_paths.profile_db) $(db_paths.seq_db) $(db_paths.realigned_result_db) $hits_tsv --format-output query,target,tseq`,
            "convertalis", logs_dir)

        all_hits, filtered_hits = _collect_hits(hits_tsv, seed_set)
        raw_hits_fasta = joinpath(run_dir, "mmseqs_hits_raw.fasta")
        filtered_fasta = joinpath(hmm_dir, "mmseqs_hits_filtered.fasta")
        write_fasta(raw_hits_fasta, all_hits)
        write_fasta(filtered_fasta, filtered_hits)

        hmm_path = joinpath(seed_dir, "seed.hmm")
        annotated_seed = joinpath(seed_dir, "seed_annotated.sto")
        0.0 <= hmmbuild_symfrac <= 1.0 ||
            error("hmmbuild_symfrac must be between 0.0 and 1.0.")
        _run_labeled(
            `$(HMMER_jll.hmmbuild()) --symfrac $(Float64(hmmbuild_symfrac)) -O $annotated_seed $hmm_path $sanitized_seed`,
            "hmmbuild", logs_dir)
        normalize_stockholm_annotations!(annotated_seed)

        aligned_sto = joinpath(hmm_dir, "alignment_with_hits.sto")
        hits_aligned_sto = joinpath(hmm_dir, "aligned_hits_only.sto")
        if isempty(filtered_hits)
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

        match_stockholm = joinpath(unpack_dir, "$(target.transcript_id)_matchonly.sto")
        open(match_stockholm, "w") do io
            run(pipeline(`$(HMMER_jll.esl_alimask()) --rf-is-mask $aligned_sto`, stdout = io))
        end
        normalize_stockholm_annotations!(match_stockholm)

        full_stockholm = joinpath(unpack_dir, "$(target.transcript_id)_full.sto")
        full_alignment = _reorder_alignment(read_file(aligned_sto, Stockholm; keepinserts = true), seed_names)
        match_alignment = _reorder_alignment(
            read_file(match_stockholm, Stockholm; keepinserts = true), seed_names)
        write_file(full_stockholm, full_alignment, Stockholm)
        write_file(match_stockholm, match_alignment, Stockholm)

        a3m_path = joinpath(unpack_dir, "$(target.transcript_id)_expanded.a3m")
        write_file(a3m_path, match_alignment, A3M)
        hits_copy = joinpath(unpack_dir, "$(target.transcript_id)_hits_raw.fasta")
        cp(raw_hits_fasta, hits_copy; force = true)

        return ExpansionResult(;
            run_dir,
            seed_stockholm = archived_seed_sto,
            seed_fasta = archived_seed_fasta,
            hits_fasta = hits_copy,
            full_stockholm,
            match_stockholm,
            a3m_path,
            db_dir,
            hmm_dir,
            logs_dir,
            n_hits = length(all_hits),
            n_new_hits = length(filtered_hits),
            status = :ok
        )
    finally
        isdir(tmp_dir) && rm(tmp_dir; recursive = true, force = true)
    end
end

end
