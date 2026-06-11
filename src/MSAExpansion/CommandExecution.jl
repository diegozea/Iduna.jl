# Command Execution
# -----------------

function _run_labeled(cmd::Cmd, label::AbstractString, logs_dir::AbstractString)
    mkpath(logs_dir)
    @debug "Running MSA expansion command." label
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
        mmseqs_threads::Union{Nothing, Integer} = nothing)
    mmseqs_bin = MMseqs2_jll.mmseqs()
    mkpath(db_dir)
    seed_db = joinpath(db_dir, "seed_msa_db")
    profile_db = joinpath(db_dir, "seed_profile_db")
    search_result_db = joinpath(db_dir, "search_results_db")
    expanded_result_db = joinpath(db_dir, "expanded_results_db")
    realigned_result_db = joinpath(db_dir, "realigned_results_db")
    aln_db = string(base_db, "_aln")
    seq_db = string(base_db, "_seq")

    # Build a profile from the seed MSA and search it against the MMseqs DB.
    _run_labeled(`$(mmseqs_bin) convertmsa $seed_msa_sto $seed_db`, "convertmsa", logs_dir)

    profile_cmd = `$(mmseqs_bin) msa2profile $seed_db $profile_db --match-mode $(match_mode)`
    match_ratio !== nothing &&
        (profile_cmd = `$profile_cmd --match-ratio $(Float64(match_ratio))`)
    _run_labeled(profile_cmd, "msa2profile", logs_dir)

    search_cmd = `$(mmseqs_bin) search $profile_db $base_db $search_result_db $tmp_dir -a`
    mmseqs_threads !== nothing &&
        (search_cmd = `$search_cmd --threads $(Int(mmseqs_threads))`)
    _run_labeled(search_cmd, "search", logs_dir)

    # The search results are centroid/consensus hits. `expandaln` expands those
    # hits to all cluster members, which remains Iduna's main MSA output.
    expandaln_cmd = `$(mmseqs_bin) expandaln $profile_db $base_db $search_result_db $aln_db $expanded_result_db`
    mmseqs_threads !== nothing &&
        (expandaln_cmd = `$expandaln_cmd --threads $(Int(mmseqs_threads))`)
    _run_labeled(expandaln_cmd, "expandaln", logs_dir)

    align_cmd = `$(mmseqs_bin) align $profile_db $seq_db $expanded_result_db $realigned_result_db`
    mmseqs_threads !== nothing &&
        (align_cmd = `$align_cmd --threads $(Int(mmseqs_threads))`)
    _run_labeled(align_cmd, "align", logs_dir)

    return (; seed_db, profile_db, seq_db, search_result_db,
        expanded_result_db, realigned_result_db)
end
