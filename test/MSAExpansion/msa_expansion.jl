@testset "MSAExpansion" begin
    mktempdir() do tmp
        source = joinpath(tmp, "in.sto")
        dest = joinpath(tmp, "out.sto")
        write(source, "seq1 AC\n")
        Iduna.MSAExpansion.prepare_stockholm_for_mmseqs(source, dest)
        text = read(dest, String)
        @test startswith(text, "# STOCKHOLM 1.0")
        @test endswith(chomp(text), "//")

        complete = joinpath(tmp, "complete.sto")
        complete_out = joinpath(tmp, "complete_out.sto")
        write(complete, "# STOCKHOLM 1.0\nseq1 AC\n//\n")
        Iduna.MSAExpansion.prepare_stockholm_for_mmseqs(complete, complete_out)
        complete_text = read(complete_out, String)
        @test complete_text == "# STOCKHOLM 1.0\nseq1 AC\n//\n"

        annotated_seed = joinpath(tmp, "annotated_seed.sto")
        sanitized_seed = joinpath(tmp, "sanitized_seed.sto")
        write(annotated_seed, """
        # STOCKHOLM 1.0
        #=GF SExonCodeMap "a"=>"1_0"
        seq1 AC
        #=GC SExonCode aa
        //
        """)
        Iduna.MSAExpansion.prepare_stockholm_for_mmseqs(annotated_seed, sanitized_seed)
        sanitized_text = read(sanitized_seed, String)
        @test !occursin("SExonCode", sanitized_text)
        @test occursin("seq1 AC", sanitized_text)

        empty_sto = joinpath(tmp, "empty.sto")
        touch(empty_sto)
        @test Iduna.MSAExpansion.normalize_stockholm_annotations!(empty_sto) == empty_sto
        @test isempty(read(empty_sto, String))

        fragmented = joinpath(tmp, "fragmented.sto")
        write(fragmented, """
        # STOCKHOLM 1.0
        # source comment
        #=GF ID family1
        #=GF DE first description
        #=GS seq1 AC accession1
        seq1 AC
        #=GR seq1 PP 99
        #=GC RF xx
        #=GF
        #=GS seq1
        seq1 DE
        #=GR seq1 PP **
        #=GC RF yy
        //
        """)
        Iduna.MSAExpansion.normalize_stockholm_annotations!(fragmented)
        normalized = read(fragmented, String)
        @test occursin("# source comment", normalized)
        @test occursin("#=GF ID family1", normalized)
        @test occursin("#=GF DE first description", normalized)
        @test occursin("#=GS seq1 AC accession1", normalized)
        @test occursin("seq1\tACDE", normalized)
        @test occursin("#=GC RF xxyy", normalized)
        @test occursin("#=GR seq1 PP 99**", normalized)

        hits_tsv = joinpath(tmp, "hits.tsv")
        write(hits_tsv, "query\tseed one\tACD-\nquery\thit one\tACDF\n")
        all_hits, filtered_hits = Iduna.MSAExpansion._collect_hits(hits_tsv, Set(["seed"]))
        @test all_hits == [("seed", "ACD"), ("hit", "ACDF")]
        @test filtered_hits == [("hit", "ACDF")]

        noisy_hits_tsv = joinpath(tmp, "noisy_hits.tsv")
        write(noisy_hits_tsv, """
        too-short
        query\t\tACDE
        query\thit one\t---
        query\thit one\tACD-
        query\thit one\tACDE
        query\tseed one\tACDF
        query\tnew_hit\tacdg
        """)
        noisy_all_hits,
        noisy_filtered_hits = Iduna.MSAExpansion._collect_hits(noisy_hits_tsv, Set(["seed"]))
        @test noisy_all_hits ==
              [("hit", "ACD"), ("seed", "ACDF"), ("new_hit", "ACDG")]
        @test noisy_filtered_hits == [("hit", "ACD"), ("new_hit", "ACDG")]

        logs_dir = joinpath(tmp, "run_logs")
        @test Iduna.MSAExpansion._run_labeled(
            `sh -c "printf msa-out; printf msa-err >&2"`, "mock", logs_dir) === nothing
        @test read(joinpath(logs_dir, "mock_stdout.log"), String) == "msa-out"
        @test read(joinpath(logs_dir, "mock_stderr.log"), String) == "msa-err"

        empty_hits = joinpath(tmp, "empty_hits.fasta")
        touch(empty_hits)
        @test Iduna.MSAExpansion._cached_hit_counts(empty_hits, Set(["seed"])) ==
              (n_hits = 0, n_new_hits = 0)

        reorder_sto = joinpath(tmp, "reorder.sto")
        write(reorder_sto, "# STOCKHOLM 1.0\nb AC\na AC\nc AC\n//\n")
        reordered = Iduna.MSAExpansion._reorder_alignment(
            Iduna.MSAExpansion.read_file(reorder_sto, Iduna.MSAExpansion.Stockholm),
            ["a", "missing"])
        @test String.(Iduna.MSAExpansion.sequencenames(reordered)) == ["a", "b", "c"]

        sexon_sto = joinpath(tmp, "sexon_project.sto")
        write(sexon_sto, """
        # STOCKHOLM 1.0
        seed ACdeFG
        hit AC-eFG
        #=GC RF xxxxxx
        //
        """)
        sexon_msa = Iduna.MSAExpansion.read_file(
            sexon_sto, Iduna.MSAExpansion.Stockholm; keepinserts = true)
        @test Iduna.MSAExpansion.getannotcolumn(sexon_msa, "Aligned", "") == "110011"
        archived = (;
            seed_s_exon_codes = "0β23",
            seed_match_s_exon_codes = "0β23",
            seed_s_exon_code_map = [
                '0' => "1_0", 'β' => "2_0", '2' => "3_0", '3' => "4_0"]
        )
        Iduna.MSAExpansion._restore_s_exon_annotations!(sexon_msa, archived)
        @test Iduna.Utils.s_exon_codes(sexon_msa) == "0β..23"
        @test Iduna.Utils.s_exon_code_map(sexon_msa)['3'] == "4_0"

        annotated_seed_sto = joinpath(tmp, "annotated_seed_rf.sto")
        write(annotated_seed_sto, """
        # STOCKHOLM 1.0
        seed ACdE
        #=GC RF xxxx
        //
        """)
        annotated_seed = Iduna.MSAExpansion.read_file(
            annotated_seed_sto, Iduna.MSAExpansion.Stockholm; keepinserts = true)
        @test Iduna.MSAExpansion.getannotcolumn(annotated_seed, "Aligned", "") == "1101"
        @test Iduna.MSAExpansion._seed_match_s_exon_codes("0123", annotated_seed) ==
              "013"
        masked_archived = (;
            seed_s_exon_codes = "0123",
            seed_match_s_exon_codes = "013",
            seed_s_exon_code_map = ['0' => "1_0", '1' => "2_0", '3' => "4_0"]
        )
        Iduna.MSAExpansion._restore_s_exon_annotations!(sexon_msa, masked_archived)
        @test Iduna.Utils.s_exon_codes(sexon_msa) == "01..3."

        @test Iduna.MSAExpansion._rf_match_state_mask("xx..xx", 6) ==
              [true, true, false, false, true, true]
        @test Iduna.MSAExpansion._rf_match_state_mask("", 6) === nothing
        @test Iduna.MSAExpansion._aligned_match_state_mask("", 6) === nothing
        plain_sto = joinpath(tmp, "plain_no_annotations.sto")
        write(plain_sto, "# STOCKHOLM 1.0\nseed AC\n//\n")
        plain_msa = Iduna.MSAExpansion.read_file(
            plain_sto, Iduna.MSAExpansion.Stockholm; keepinserts = true)
        @test Iduna.MSAExpansion._match_state_mask(
            plain_msa; default_aligned = true) == [true, true]
        plain_msa_without_insert_annotation = Iduna.MSAExpansion.read_file(
            plain_sto, Iduna.MSAExpansion.Stockholm)
        @test Iduna.MSAExpansion._match_state_mask(
            plain_msa_without_insert_annotation; default_aligned = false) ==
              [false, false]

        gene_id = "ENSG00000198821"
        transcript_id = "ENST00000362089.10"
        seed_sto = joinpath(tmp, "seed.sto")
        write(seed_sto, "# STOCKHOLM 1.0\nseed ACDE\n//\n")
        seed = Iduna.SeedSelection(;
            pid = 10.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = seed_sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        target = Iduna.ResolvedTarget(;
            input_id = transcript_id,
            input_kind = :ensembl_transcript,
            ensembl_gene_id = gene_id,
            transcript_id
        )

        db = joinpath(tmp, "mock_mmseqs_db")
        touch("$(db).dbtype")
        touch("$(db)_aln.dbtype")
        touch("$(db)_seq.dbtype")

        @test_throws ErrorException Iduna.MSAExpansion.expand_msa(
            target, seed, joinpath(tmp, "invalid_symfrac"); mmseqs_db = db,
            hmmbuild_symfrac = 1.1)
        @test !isdir(joinpath(tmp, "invalid_symfrac", "expansion"))

        function write_outputs(outputs; centroids::Bool = false)
            mkpath(dirname(outputs.full_stockholm))
            write(outputs.full_stockholm,
                "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
            write(outputs.match_stockholm,
                "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
            write(outputs.a3m_path, "cached a3m\n")
            write(outputs.hits_fasta, ">seed one\nACDE\n>hit one\nACDF\n>hit_two\nACDG\n")
            if centroids
                mkpath(dirname(outputs.centroid_full_stockholm))
                write(outputs.centroid_full_stockholm,
                    "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
                write(outputs.centroid_match_stockholm,
                    "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
                write(outputs.centroid_a3m_path, "cached centroid a3m\n")
                write(outputs.centroid_hits_fasta, ">hit one\nACDF\n")
            end
            return outputs
        end

        function build_self_hit_mmseqs_db(root::AbstractString)
            mmseqs = Iduna.MSAExpansion.MMseqs2_jll.mmseqs()
            mkpath(root)
            fasta = joinpath(root, "self_hit.fasta")
            write(fasta, ">seed\nACDEFGHIKLMNPQRSTVWY\n")
            db = joinpath(root, "self_hit_db")
            seq_db = string(db, "_seq")
            aln_db = string(db, "_aln")
            tmp_dir = joinpath(root, "tmp")
            logs_dir = joinpath(root, "logs")
            mkpath.((tmp_dir, logs_dir))
            Iduna.MSAExpansion._run_labeled(
                `$(mmseqs) createdb $fasta $db`, "createdb", logs_dir)
            Iduna.MSAExpansion._run_labeled(
                `$(mmseqs) createdb $fasta $seq_db`, "createdb_seq", logs_dir)
            Iduna.MSAExpansion._run_labeled(
                `$(mmseqs) search $db $seq_db $aln_db $tmp_dir -a --threads 1`,
                "search", logs_dir)
            return db
        end

        archive_seed_sto = joinpath(tmp, "archive_seed.sto")
        archive_seed_fasta = joinpath(tmp, "archive_seed.fasta")
        write(archive_seed_sto,
            "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed AC\n#=GC SExonCode 00\n//\n")
        write(archive_seed_fasta, ">seed\nAC\n")
        archive_ctx = (;
            seed_dir = joinpath(tmp, "archive_seeds"),
            seed_stockholm = archive_seed_sto,
            seed_fasta = archive_seed_fasta)
        mkpath(archive_ctx.seed_dir)
        archived_with_fasta = Iduna.MSAExpansion._archive_expansion_seed(
            seed, archive_ctx)
        @test isfile(archived_with_fasta.archived_seed_fasta)
        @test read(archived_with_fasta.archived_seed_fasta, String) == ">seed\nAC\n"
        @test archived_with_fasta.seed_s_exon_codes == "00"

        run_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_10.00")
        outputs = write_outputs(
            Iduna.MSAExpansion._expansion_output_paths(run_dir, transcript_id))
        identity = Iduna.MSAExpansion._expansion_identity(
            target, seed, seed_sto, nothing, db;
            match_mode = 1,
            match_ratio = nothing,
            hmmbuild_symfrac = 0.0,
            centroids = false)
        Iduna.MSAExpansion._write_step_state(run_dir, :done, identity, outputs)
        state = Iduna.MSAExpansion._read_step_state(run_dir)
        @test state.status == "done"
        @test state.step == "msa_expansion"
        @test Iduna.MSAExpansion._step_state_unreadable_message(nothing) ==
              "state file disappeared while reading"
        hits_fasta = outputs.hits_fasta

        cached = @test_logs (:info, r"Reusing cached MSA expansion") match_mode=:any Iduna.MSAExpansion.expand_msa(
            target, seed, tmp; mmseqs_db = db, threads = 1)
        cached_with_new_threads = Iduna.MSAExpansion.expand_msa(
            target, seed, tmp; mmseqs_db = db, threads = 8)
        hits_msa = Iduna.MSAExpansion.read_file(hits_fasta, Iduna.MSAExpansion.FASTA)
        @test cached.status === :skipped
        @test cached_with_new_threads.status === :skipped
        @test cached.n_hits == Iduna.MSAExpansion.nsequences(hits_msa)
        @test cached.n_hits == 3
        @test cached.n_new_hits == 2

        changed_mode_identity = Iduna.MSAExpansion._expansion_identity(
            target, seed, seed_sto, nothing, db;
            match_mode = 3,
            match_ratio = nothing,
            hmmbuild_symfrac = 0.0,
            centroids = false)
        changed_mode = Iduna.MSAExpansion._classify_step_state(
            run_dir, changed_mode_identity, outputs)
        @test changed_mode.reusable === false
        @test changed_mode.status === :outdated

        changed_ratio_identity = Iduna.MSAExpansion._expansion_identity(
            target, seed, seed_sto, nothing, db;
            match_mode = 1,
            match_ratio = 0.75,
            hmmbuild_symfrac = 0.0,
            centroids = false)
        changed_ratio = Iduna.MSAExpansion._classify_step_state(
            run_dir, changed_ratio_identity, outputs)
        @test changed_ratio.reusable === false
        @test changed_ratio.status === :outdated

        changed_symfrac_identity = Iduna.MSAExpansion._expansion_identity(
            target, seed, seed_sto, nothing, db;
            match_mode = 1,
            match_ratio = nothing,
            hmmbuild_symfrac = 0.5,
            centroids = false)
        changed_symfrac = Iduna.MSAExpansion._classify_step_state(
            run_dir, changed_symfrac_identity, outputs)
        @test changed_symfrac.reusable === false
        @test changed_symfrac.status === :outdated

        write(seed_sto, "# STOCKHOLM 1.0\nseed ACDF\n//\n")
        changed_seed_identity = Iduna.MSAExpansion._expansion_identity(
            target, seed, seed_sto, nothing, db;
            match_mode = 1,
            match_ratio = nothing,
            hmmbuild_symfrac = 0.0,
            centroids = false)
        changed_seed = Iduna.MSAExpansion._classify_step_state(
            run_dir, changed_seed_identity, outputs)
        @test changed_seed.reusable === false
        @test changed_seed.status === :outdated
        write(seed_sto, "# STOCKHOLM 1.0\nseed ACDE\n//\n")

        centroids_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_20.00")
        centroids_seed = Iduna.SeedSelection(;
            pid = 20.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = seed_sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        centroids_outputs = write_outputs(
            Iduna.MSAExpansion._expansion_output_paths(
                centroids_dir, transcript_id; centroids = true);
            centroids = true)
        centroids_identity = Iduna.MSAExpansion._expansion_identity(
            target, centroids_seed, seed_sto, nothing, db;
            match_mode = 1,
            match_ratio = nothing,
            hmmbuild_symfrac = 0.0,
            centroids = true)
        Iduna.MSAExpansion._write_step_state(
            centroids_dir, :done, centroids_identity, centroids_outputs)
        write(hits_fasta, ">seed one\nACDE\n>hit one\nACDF\n>hit_two\nACDG\n")

        cached_centroids = Iduna.MSAExpansion.expand_msa(
            target, centroids_seed, tmp; mmseqs_db = db, centroids = true)
        @test cached_centroids.status === :skipped
        @test isfile(centroids_outputs.s_exon_blocks_tsv)
        @test isfile(centroids_outputs.centroid_s_exon_blocks_tsv)

        missing_centroid_dir = joinpath(
            tmp, "expansion", gene_id, transcript_id, "pid_30.00")
        missing_centroid_outputs = write_outputs(
            Iduna.MSAExpansion._expansion_output_paths(
            missing_centroid_dir, transcript_id; centroids = true))
        missing_centroid_seed = Iduna.SeedSelection(;
            pid = 30.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = seed_sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        missing_centroid_identity = Iduna.MSAExpansion._expansion_identity(
            target, missing_centroid_seed, seed_sto, nothing, db;
            match_mode = 1,
            match_ratio = nothing,
            hmmbuild_symfrac = 0.0,
            centroids = true)
        Iduna.MSAExpansion._write_step_state(
            missing_centroid_dir, :done, missing_centroid_identity,
            missing_centroid_outputs)
        missing_centroids = Iduna.MSAExpansion._classify_step_state(
            missing_centroid_dir, missing_centroid_identity, missing_centroid_outputs)
        @test missing_centroids.reusable === false
        @test missing_centroids.status === :unfinished

        incomplete_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_35.00")
        incomplete_outputs = Iduna.MSAExpansion._expansion_output_paths(
            incomplete_dir, transcript_id)
        mkpath(incomplete_dir)
        incomplete = Iduna.MSAExpansion._classify_step_state(
            incomplete_dir, identity, incomplete_outputs)
        @test incomplete.reusable === false
        @test incomplete.status === :unfinished
        @test occursin("incomplete outputs", incomplete.warning)

        missing_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_36.00")
        missing_outputs = Iduna.MSAExpansion._expansion_output_paths(
            missing_dir, transcript_id)
        missing = Iduna.MSAExpansion._classify_step_state(
            missing_dir, identity, missing_outputs)
        @test missing.reusable === false
        @test missing.status === :missing
        @test missing.warning === nothing

        unreadable_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_37.00")
        unreadable_outputs = Iduna.MSAExpansion._expansion_output_paths(
            unreadable_dir, transcript_id)
        mkpath(unreadable_dir)
        write(Iduna.MSAExpansion._step_state_path(unreadable_dir), "{not json")
        unreadable = Iduna.MSAExpansion._classify_step_state(
            unreadable_dir, identity, unreadable_outputs)
        @test unreadable.reusable === false
        @test unreadable.status === :outdated
        @test occursin("Could not read MSA expansion step_state.json", unreadable.warning)

        legacy_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_40.00")
        legacy_outputs = write_outputs(
            Iduna.MSAExpansion._expansion_output_paths(legacy_dir, transcript_id))
        legacy_seed = Iduna.SeedSelection(;
            pid = 40.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = seed_sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        legacy_identity = Iduna.MSAExpansion._expansion_identity(
            target, legacy_seed, seed_sto, nothing, db;
            match_mode = 1,
            match_ratio = nothing,
            hmmbuild_symfrac = 0.0,
            centroids = false)
        legacy = Iduna.MSAExpansion._classify_step_state(
            legacy_dir, legacy_identity, legacy_outputs)
        @test legacy.reusable === false
        @test legacy.status === :outdated
        @test occursin("no step_state.json", legacy.warning)

        stale_failure_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_50.00")
        stale_failure_outputs = write_outputs(
            Iduna.MSAExpansion._expansion_output_paths(stale_failure_dir, transcript_id))
        stale_failure_seed = Iduna.SeedSelection(;
            pid = 50.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = seed_sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        @test_throws Base.ProcessFailedException Iduna.MSAExpansion.expand_msa(
            target, stale_failure_seed, tmp; mmseqs_db = db, threads = 1)
        stale_failure_state = Iduna.MSAExpansion._read_step_state(stale_failure_dir)
        @test stale_failure_state.status == "failed"
        @test occursin("no step_state.json", stale_failure_state.warnings[1])
        @test occursin("ProcessFailedException", stale_failure_state.exception.type)
        cache_warning = read(joinpath(stale_failure_dir, "logs", "cache_warning.log"),
            String)
        @test occursin("no step_state.json", cache_warning)
        @test !isfile(stale_failure_outputs.full_stockholm)

        success_seed_sto = joinpath(tmp, "success_seed.sto")
        write(success_seed_sto, "# STOCKHOLM 1.0\nseed ACDEFGHIKLMNPQRSTVWY\n//\n")
        success_seed = Iduna.SeedSelection(;
            pid = 60.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = success_seed_sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        success_db = build_self_hit_mmseqs_db(joinpath(tmp, "self_hit_mmseqs"))
        success = Iduna.MSAExpansion.expand_msa(
            target, success_seed, tmp; mmseqs_db = success_db, centroids = true,
            threads = 1)
        @test success.status === :ok
        @test success.n_new_hits == 0
        @test isfile(joinpath(success.run_dir, "centroid_msa",
            "$(transcript_id)_centroids.a3m"))
        success_state = Iduna.MSAExpansion._read_step_state(success.run_dir)
        @test success_state.status == "done"

        Iduna.MSAExpansion._write_step_state(legacy_dir, :failed, legacy_identity,
            legacy_outputs; exception = (; type = "ErrorException", message = "failed"))
        failed = Iduna.MSAExpansion._classify_step_state(
            legacy_dir, legacy_identity, legacy_outputs)
        @test failed.reusable === false
        @test failed.status === :failed

        Iduna.MSAExpansion._write_step_state(legacy_dir, :unfinished, legacy_identity,
            legacy_outputs)
        unfinished = Iduna.MSAExpansion._classify_step_state(
            legacy_dir, legacy_identity, legacy_outputs)
        @test unfinished.reusable === false
        @test unfinished.status === :unfinished
    end
end
