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
            write(outputs.full_stockholm, "cached full\n")
            write(outputs.match_stockholm, "cached match\n")
            write(outputs.a3m_path, "cached a3m\n")
            write(outputs.hits_fasta, ">seed one\nACDE\n>hit one\nACDF\n>hit_two\nACDG\n")
            if centroids
                mkpath(dirname(outputs.centroid_full_stockholm))
                write(outputs.centroid_full_stockholm, "cached centroid full\n")
                write(outputs.centroid_match_stockholm, "cached centroid match\n")
                write(outputs.centroid_a3m_path, "cached centroid a3m\n")
                write(outputs.centroid_hits_fasta, ">hit one\nACDF\n")
            end
            return outputs
        end

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
        hits_fasta = outputs.hits_fasta

        cached = Iduna.MSAExpansion.expand_msa(target, seed, tmp; mmseqs_db = db,
            threads = 1)
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
