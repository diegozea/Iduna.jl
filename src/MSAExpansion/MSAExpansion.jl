module MSAExpansion

import HMMER_jll
import MMseqs2_jll

using Dates: UTC, now
using MIToS.MSA: A3M, AbstractMultipleSequenceAlignment, FASTASequences, Stockholm,
                 getannotcolumn, ncolumns, nsequences, read_file, sequence_id,
                 sequencenames, write_file

using ..Utils: ExpansionResult, ResolvedTarget, SeedSelection, _resolve_artifact_path,
               ensure_mmseqs_db, format_pid, format_pid_dir, run_logged, safe_rm,
               s_exon_code_map, s_exon_codes, set_s_exon_annotations!, write_fasta,
               write_s_exon_blocks_tsv, _classify_stage_state, _file_sha256,
               _read_stage_state, _stage_state_path, _write_stage_state

export expand_msa,
       normalize_stockholm_annotations!,
       prepare_stockholm_for_mmseqs

_normalize_id(s::AbstractString) = String(split(String(s))[1])

const _STEP_STATE_FILE = "stage_state.json"

function _expansion_root(workdir::AbstractString)
    joinpath(workdir, "expansion")
end

function _expansion_output_paths(run_dir::AbstractString,
        transcript_id::AbstractString;
        centroids::Bool = false)
    expanded_dir = joinpath(run_dir, "expanded_msa")
    outputs = (;
        full_stockholm = joinpath(expanded_dir, "$(transcript_id)_full.sto"),
        match_stockholm = joinpath(expanded_dir, "$(transcript_id)_matchonly.sto"),
        a3m_path = joinpath(expanded_dir, "$(transcript_id)_expanded.a3m"),
        s_exon_blocks_tsv = joinpath(expanded_dir, "$(transcript_id)_s_exon_blocks.tsv"),
        hits_fasta = joinpath(expanded_dir, "$(transcript_id)_hits_raw.fasta")
    )
    if !centroids
        return outputs
    end
    centroid_dir = joinpath(run_dir, "centroid_msa")
    return merge(outputs,
        (;
            centroid_full_stockholm = joinpath(
                centroid_dir, "$(transcript_id)_centroids_full.sto"),
            centroid_match_stockholm = joinpath(
                centroid_dir, "$(transcript_id)_centroids_matchonly.sto"),
            centroid_a3m_path = joinpath(centroid_dir, "$(transcript_id)_centroids.a3m"),
            centroid_s_exon_blocks_tsv = joinpath(
                centroid_dir, "$(transcript_id)_centroids_s_exon_blocks.tsv"),
            centroid_hits_fasta = joinpath(
                centroid_dir, "$(transcript_id)_centroid_hits_raw.fasta")
        ))
end

_step_state_path(run_dir::AbstractString) = _stage_state_path(run_dir)

function _required_expansion_outputs(outputs::NamedTuple)
    return Dict(String(name) => path
    for (name, path) in pairs(outputs)
    if !endswith(String(name), "s_exon_blocks_tsv"))
end

function _expansion_identity(target::ResolvedTarget,
        seed::SeedSelection,
        seed_stockholm::AbstractString,
        seed_fasta::Union{Nothing, AbstractString},
        mmseqs_db::AbstractString;
        match_mode::Integer,
        match_ratio::Union{Nothing, Real},
        hmmbuild_symfrac::Real,
        centroids::Bool)
    return (;
        target = (;
            input_id = target.input_id,
            input_kind = String(target.input_kind),
            uniprot_id = target.uniprot_id,
            ensembl_gene_id = target.ensembl_gene_id,
            transcript_id = target.transcript_id,
            ensembl_protein_id = target.ensembl_protein_id,
            species = target.species
        ),
        seed = (;
            pid = Float64(seed.pid),
            stockholm_sha256 = _file_sha256(seed_stockholm),
            fasta_sha256 = seed_fasta === nothing || !isfile(seed_fasta) ? nothing :
                           _file_sha256(seed_fasta)
        ),
        expansion = (;
            mmseqs_db = abspath(String(mmseqs_db)),
            match_mode = Int(match_mode),
            match_ratio = match_ratio === nothing ? nothing : Float64(match_ratio),
            hmmbuild_symfrac = Float64(hmmbuild_symfrac),
            centroids = Bool(centroids)
        )
    )
end

function _exception_summary(err)
    return (;
        type = string(typeof(err)),
        message = sprint(showerror, err)
    )
end

function _expansion_stage_key(identity)
    return "expansion:$(identity.target.ensembl_gene_id):$(identity.target.transcript_id):$(format_pid_dir(identity.seed.pid))"
end

function _write_step_state(run_dir::AbstractString,
        status::Symbol,
        identity,
        outputs::NamedTuple;
        warnings::AbstractVector{<:AbstractString} = String[],
        exception = nothing,
        action::Union{Nothing, Symbol, AbstractString} = nothing)
    return _write_stage_state(run_dir;
        stage = "msa_expansion",
        stage_key = _expansion_stage_key(identity),
        status,
        identity,
        outputs,
        warnings,
        exception,
        action,
        extra = (; step = "msa_expansion"))
end

_read_step_state(run_dir::AbstractString) = _read_stage_state(run_dir)

function _step_state_unreadable_message(state)
    if state === nothing
        return "state file disappeared while reading"
    elseif state isa NamedTuple && haskey(state, :unreadable)
        return state.unreadable
    end
    return nothing
end

function _classify_step_state(run_dir::AbstractString, identity, outputs::NamedTuple)
    return _classify_stage_state(run_dir, identity, _required_expansion_outputs(outputs);
        stage_label = "MSA expansion")
end

function _write_cache_warning(logs_dir::AbstractString, warning::AbstractString)
    mkpath(logs_dir)
    open(joinpath(logs_dir, "cache_warning.log"), "a") do io
        println(io, "[", now(UTC), "] ", warning)
    end
    return nothing
end

function _is_s_exon_provenance_stockholm_line(line::AbstractString)
    return startswith(line, "#=GC SExonCode") ||
           startswith(line, "#=GF SExonCodeMap")
end

function prepare_stockholm_for_mmseqs(source::AbstractString, dest::AbstractString)
    mkpath(dirname(dest))
    lines = readlines(source)
    has_header = !isempty(lines) && startswith(strip(lines[1]), "# STOCKHOLM")
    open(dest, "w") do io
        # MMseqs expects a complete Stockholm file with header and terminator.
        has_header || println(io, "# STOCKHOLM 1.0")
        for line in lines
            _is_s_exon_provenance_stockholm_line(line) && continue
            println(io, line)
        end
        if isempty(lines) || strip(lines[end]) != "//"
            println(io, "//")
        end
    end
    return dest
end

function _stockholm_annotation_data()
    return (;
        comments = String[],
        gf_order = String[],
        gf_data = Dict{String, Vector{String}}(),
        gs_order = Tuple{String, String}[],
        gs_data = Dict{Tuple{String, String}, Vector{String}}(),
        gc_order = String[],
        gc_data = Dict{String, String}(),
        gr_order = Tuple{String, String}[],
        gr_data = Dict{Tuple{String, String}, String}(),
        seq_order = String[],
        seq_data = Dict{String, String}()
    )
end

function _append_stockholm_vector!(order, data, key, value)
    haskey(data, key) || (push!(order, key); data[key] = String[])
    push!(data[key], value)
    return nothing
end

function _append_stockholm_string!(order, data, key, value)
    if haskey(data, key)
        data[key] = string(data[key], value)
    else
        push!(order, key)
        data[key] = value
    end
    return nothing
end

function _record_gf_line!(records, line::AbstractString)
    parts = split(line; limit = 3)
    length(parts) < 3 && return nothing
    _append_stockholm_vector!(records.gf_order, records.gf_data, parts[2], parts[3])
    return nothing
end

function _record_gs_line!(records, line::AbstractString)
    parts = split(line; limit = 4)
    length(parts) < 4 && return nothing
    key = (parts[2], parts[3])
    _append_stockholm_vector!(records.gs_order, records.gs_data, key, parts[4])
    return nothing
end

function _record_gc_line!(records, line::AbstractString)
    parts = split(line; limit = 3)
    length(parts) < 3 && return nothing
    _append_stockholm_string!(
        records.gc_order, records.gc_data, parts[2], replace(parts[3], ' ' => ""))
    return nothing
end

function _record_gr_line!(records, line::AbstractString)
    parts = split(line; limit = 4)
    length(parts) < 4 && return nothing
    key = (parts[2], parts[3])
    _append_stockholm_string!(
        records.gr_order, records.gr_data, key, replace(parts[4], ' ' => ""))
    return nothing
end

function _record_sequence_line!(records, line::AbstractString)
    parts = split(line; limit = 2)
    length(parts) < 2 && return nothing
    name = strip(parts[1])
    fragment = replace(strip(parts[2]), ' ' => "")
    _append_stockholm_string!(records.seq_order, records.seq_data, name, fragment)
    return nothing
end

function _record_stockholm_line!(records, line::AbstractString)
    stripped = strip(line)
    isempty(stripped) && return nothing
    (startswith(stripped, "# STOCKHOLM") || startswith(stripped, "//")) &&
        return nothing
    startswith(line, "#=GF") && return _record_gf_line!(records, line)
    startswith(line, "#=GS") && return _record_gs_line!(records, line)
    startswith(line, "#=GC") && return _record_gc_line!(records, line)
    startswith(line, "#=GR") && return _record_gr_line!(records, line)
    startswith(line, '#') && (push!(records.comments, stripped); return nothing)
    return _record_sequence_line!(records, line)
end

function _write_stockholm_metadata(io, records)
    foreach(line -> println(io, line), records.comments)
    for feature in records.gf_order, value in records.gf_data[feature]

        println(io, "#=GF ", feature, ' ', value)
    end
    for key in records.gs_order, value in records.gs_data[key]

        println(io, "#=GS ", key[1], ' ', key[2], ' ', value)
    end
    for feature in records.gc_order
        println(io, "#=GC ", feature, ' ', records.gc_data[feature])
    end
    return nothing
end

function _write_stockholm_sequences(io, records)
    for name in records.seq_order
        println(io, name, '\t', records.seq_data[name])
        for key in records.gr_order
            key[1] == name || continue
            println(io, "#=GR ", key[1], ' ', key[2], ' ', records.gr_data[key])
        end
    end
    return nothing
end

function _write_normalized_stockholm(path::AbstractString, records)
    open(path, "w") do io
        println(io, "# STOCKHOLM 1.0")
        _write_stockholm_metadata(io, records)
        _write_stockholm_sequences(io, records)
        println(io, "//")
    end
    return path
end

function normalize_stockholm_annotations!(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return path

    records = _stockholm_annotation_data()
    # Merge split Stockholm records while keeping the original output order.
    for line in lines
        _record_stockholm_line!(records, line)
    end
    return _write_normalized_stockholm(path, records)
end

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

function _is_rf_match(rf_char::Char)
    return !(rf_char == '.' || rf_char == '-' || rf_char == ' ')
end

function _rf_match_state_mask(rf::AbstractString, ncols::Integer)
    if isempty(rf)
        return nothing
    end
    length(rf) == ncols ||
        error("RF annotation has $(length(rf)) characters, but the MSA has $(ncols) columns.")
    return [_is_rf_match(rf_char) for rf_char in rf]
end

function _aligned_match_state_mask(aligned::AbstractString, ncols::Integer)
    if isempty(aligned)
        return nothing
    end
    length(aligned) == ncols ||
        error("Aligned annotation has $(length(aligned)) characters, but the MSA has $(ncols) columns.")
    return [aligned_char == '1' ? true :
            aligned_char == '0' ? false :
            error("Aligned annotation contains $(repr(aligned_char)); expected '0' or '1'.")
            for aligned_char in aligned]
end

function _match_state_mask(msa::AbstractMultipleSequenceAlignment;
        default_aligned::Bool)
    ncols = ncolumns(msa)
    aligned_mask = _aligned_match_state_mask(getannotcolumn(msa, "Aligned", ""), ncols)
    aligned_mask === nothing || return aligned_mask
    rf_mask = _rf_match_state_mask(getannotcolumn(msa, "RF", ""), ncols)
    rf_mask === nothing || return rf_mask
    return fill(default_aligned, ncols)
end

function _project_s_exon_codes(msa::AbstractMultipleSequenceAlignment,
        seed_match_codes::AbstractString)
    match_mask = _match_state_mask(msa; default_aligned = false)
    seed_match_chars = collect(seed_match_codes)
    io = IOBuffer(; sizehint = ncolumns(msa))
    seed_idx = 0
    for is_match in match_mask
        if is_match
            seed_idx += 1
            write(io, seed_idx <= length(seed_match_chars) ? seed_match_chars[seed_idx] :
                      '.')
        else
            write(io, '.')
        end
    end
    return String(take!(io))
end

function _seed_match_s_exon_codes(seed_codes::AbstractString,
        annotated_seed::AbstractMultipleSequenceAlignment)
    match_mask = _match_state_mask(annotated_seed; default_aligned = true)
    length(match_mask) == length(seed_codes) ||
        error("Annotated seed match mask has $(length(match_mask)) characters, but SExonCode has $(length(seed_codes)) characters.")
    io = IOBuffer(; sizehint = length(seed_codes))
    for (is_match, code) in zip(match_mask, seed_codes)
        is_match && write(io, code)
    end
    return String(take!(io))
end

function _with_seed_match_s_exon_codes(archived, annotated_seed_path::AbstractString)
    annotated_seed = read_file(annotated_seed_path, Stockholm; keepinserts = true)
    return merge(archived,
        (;
            seed_match_s_exon_codes = _seed_match_s_exon_codes(
            archived.seed_s_exon_codes, annotated_seed)))
end

function _restore_s_exon_annotations!(msa::AbstractMultipleSequenceAlignment, archived)
    seed_match_codes = :seed_match_s_exon_codes in propertynames(archived) ?
                       archived.seed_match_s_exon_codes : archived.seed_s_exon_codes
    codes = _project_s_exon_codes(msa, seed_match_codes)
    set_s_exon_annotations!(msa, codes, archived.seed_s_exon_code_map)
    return msa
end

function _write_expansion_s_exon_blocks(path::AbstractString,
        full_alignment::AbstractMultipleSequenceAlignment,
        match_alignment::AbstractMultipleSequenceAlignment,
        pid::Real)
    write_s_exon_blocks_tsv(path, match_alignment;
        alignment = "expanded_match",
        pid = Float64(pid))
    write_s_exon_blocks_tsv(path, full_alignment;
        alignment = "expanded_full",
        pid = Float64(pid),
        append = true)
    return path
end

function _ensure_alignment_s_exon_blocks(path::AbstractString,
        full_stockholm::AbstractString,
        match_stockholm::AbstractString,
        pid::Real;
        full_label::AbstractString,
        match_label::AbstractString)
    isfile(path) && return path
    (isfile(full_stockholm) && isfile(match_stockholm)) || return path
    full_alignment = read_file(full_stockholm, Stockholm; keepinserts = true)
    match_alignment = read_file(match_stockholm, Stockholm; keepinserts = true)
    write_s_exon_blocks_tsv(path, match_alignment;
        alignment = match_label,
        pid = Float64(pid))
    write_s_exon_blocks_tsv(path, full_alignment;
        alignment = full_label,
        pid = Float64(pid),
        append = true)
    return path
end

function _ensure_expansion_s_exon_blocks(outputs::NamedTuple, pid::Real)
    _ensure_alignment_s_exon_blocks(outputs.s_exon_blocks_tsv,
        outputs.full_stockholm,
        outputs.match_stockholm,
        pid;
        full_label = "expanded_full",
        match_label = "expanded_match")
    if :centroid_s_exon_blocks_tsv in propertynames(outputs)
        _ensure_alignment_s_exon_blocks(outputs.centroid_s_exon_blocks_tsv,
            outputs.centroid_full_stockholm,
            outputs.centroid_match_stockholm,
            pid;
            full_label = "centroid_full",
            match_label = "centroid_match")
    end
    return outputs.s_exon_blocks_tsv
end

function _collect_hits(hits_tsv::AbstractString, seed_set::Set{String})
    all_hits = Tuple{String, String}[]
    filtered_hits = Tuple{String, String}[]
    seen_all = Set{String}()
    seen_filtered = Set{String}()
    open(hits_tsv, "r") do io
        # Keep all hits for reporting, but align only hits not already in the seed.
        for line in eachline(io)
            parts = split(line, '\t')
            length(parts) < 3 && continue
            target_field = strip(parts[2])
            isempty(target_field) && continue
            target_name = _normalize_id(target_field)
            seq = uppercase(replace(strip(parts[3]), '-' => ""))
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
    hits = read_file(hits_fasta, FASTASequences)
    hit_names = _normalize_id.(sequence_id.(hits))
    return (;
        n_hits = length(hit_names),
        n_new_hits = count(name -> !(name in seed_set), hit_names)
    )
end

function _write_centroid_msa(transcript_id::AbstractString,
        mmseqs_db::AbstractString,
        db_paths,
        tmp_dir::AbstractString,
        centroid_dir::AbstractString,
        seed_set::Set{String},
        seed_names::Vector{String},
        sanitized_seed::AbstractString,
        hmm_path::AbstractString,
        annotated_seed::AbstractString,
        archived,
        logs_dir::AbstractString)
    mkpath(centroid_dir)
    @debug "Writing centroid MSA." transcript_id centroid_dir

    # Keep working files in the MMseqs temp area; only final centroid files go here.
    mktempdir(tmp_dir; prefix = "centroid_msa_") do centroid_tmp
        centroid_hits_tsv = joinpath(centroid_tmp, "mmseqs_centroid_hits_raw.tsv")
        # Reuse the existing search result: this saves a side MSA without
        # rerunning MMseqs search or changing the full expansion.
        _run_labeled(
            `$(MMseqs2_jll.mmseqs()) convertalis $(db_paths.profile_db) $(mmseqs_db) $(db_paths.search_result_db) $centroid_hits_tsv --format-output query,target,tseq`,
            "convertalis_centroids", logs_dir)

        # Save every centroid hit for auditing, but align only hits not already in the seed.
        all_hits, filtered_hits = _collect_hits(centroid_hits_tsv, seed_set)
        raw_hits_fasta = joinpath(centroid_dir, "$(transcript_id)_centroid_hits_raw.fasta")
        filtered_fasta = joinpath(centroid_tmp, "mmseqs_centroid_hits_filtered.fasta")
        write_fasta(raw_hits_fasta, all_hits)
        write_fasta(filtered_fasta, filtered_hits)

        # If there are no new centroid hits, the seed alignment is already the final MSA.
        aligned_sto = joinpath(centroid_tmp, "alignment_with_centroid_hits.sto")
        if isempty(filtered_hits)
            cp(annotated_seed, aligned_sto; force = true)
        else
            _run_labeled(
                `$(HMMER_jll.hmmalign()) --mapali $sanitized_seed --trim --outformat stockholm -o $aligned_sto $hmm_path $filtered_fasta`,
                "hmmalign_centroids", logs_dir)
        end
        normalize_stockholm_annotations!(aligned_sto)

        match_stockholm = joinpath(centroid_dir, "$(transcript_id)_centroids_matchonly.sto")
        open(match_stockholm, "w") do io
            run(pipeline(`$(HMMER_jll.esl_alimask()) --rf-is-mask $aligned_sto`, stdout = io))
        end
        normalize_stockholm_annotations!(match_stockholm)

        # The match-only file removes insert columns; the full file keeps them.
        full_stockholm = joinpath(centroid_dir, "$(transcript_id)_centroids_full.sto")
        full_alignment = _reorder_alignment(read_file(aligned_sto, Stockholm; keepinserts = true), seed_names)
        match_alignment = _reorder_alignment(
            read_file(match_stockholm, Stockholm; keepinserts = true), seed_names)
        _restore_s_exon_annotations!(full_alignment, archived)
        _restore_s_exon_annotations!(match_alignment, archived)
        write_file(full_stockholm, full_alignment, Stockholm)
        write_file(match_stockholm, match_alignment, Stockholm)
        write_file(joinpath(centroid_dir, "$(transcript_id)_centroids.a3m"),
            match_alignment, A3M)
        s_exon_blocks_tsv = joinpath(
            centroid_dir, "$(transcript_id)_centroids_s_exon_blocks.tsv")
        write_s_exon_blocks_tsv(s_exon_blocks_tsv, match_alignment;
            alignment = "centroid_match",
            pid = Float64(archived.seed_pid))
        write_s_exon_blocks_tsv(s_exon_blocks_tsv, full_alignment;
            alignment = "centroid_full",
            pid = Float64(archived.seed_pid),
            append = true)
    end
    return nothing
end

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
    sanitized_seed = prepare_stockholm_for_mmseqs(ctx.seed_stockholm,
        joinpath(ctx.seed_dir, "$(seed_label)_mmseqs.sto"))
    seed_alignment = read_file(archived_seed_sto, Stockholm; keepinserts = true)
    seed_names = String.(sequencenames(seed_alignment))
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

function expand_msa(target::ResolvedTarget,
        seed::SeedSelection,
        workdir::AbstractString;
        mmseqs_db::AbstractString,
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
        match_mode, match_ratio, hmmbuild_symfrac, centroids)
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

end
