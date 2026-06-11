# Hit and Centroid Outputs
# ------------------------

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
