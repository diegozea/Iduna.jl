@testset "ResultsValidation" begin
    mktempdir() do tmp
        sto = joinpath(tmp, "tiny.sto")
        write(sto, """
        # STOCKHOLM 1.0
        seq1 ACDE
        seq2 ACDF
        //
        """)
        msa = Iduna.ResultsValidation.load_msa(sto)
        @test Iduna.ResultsValidation.nsequences(msa) == 2
        @test Iduna.ResultsValidation._resolve_query_name(msa, "missing_gene", "missing_tx") ==
              "seq1"
        stats = Iduna.ResultsValidation.alignment_stats(sto)
        @test stats.n_sequences == 2
        @test stats.n_columns == 4

        a3m = joinpath(tmp, "tiny.a3m")
        write(a3m, ">seq1\nACde\n>seq2\nACfg\n")
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_msa(a3m)) == 2

        fasta = joinpath(tmp, "tiny.fasta")
        write(fasta, ">seq1\nACDE\n>seq2\nACDF\n")
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_msa(fasta)) == 2
        @test_throws ErrorException Iduna.ResultsValidation.load_msa(
            joinpath(tmp, "unknown.txt"))

        seed = Iduna.SeedSelection(;
            pid = 10.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        expansion = Iduna.ExpansionResult(;
            run_dir = tmp,
            seed_stockholm = sto,
            hits_fasta = fasta,
            full_stockholm = sto,
            match_stockholm = sto,
            a3m_path = a3m,
            db_dir = tmp,
            hmm_dir = tmp,
            logs_dir = tmp
        )
        uniprot = joinpath(tmp, "uniprot.fasta")
        write(uniprot, ">P20963\nACDF\n")
        target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001",
            uniprot_sequence_path = uniprot
        )
        validation = Iduna.ResultsValidation.validate_results(target, seed, expansion, tmp)
        @test validation.status === :warn
        @test validation.query_name == "seq1"
        @test validation.aln_mismatches == 1
        @test validation.aln_insertions == 0
        @test validation.aln_deletions == 0
        @test occursin("substitutions", only(validation.warnings))
        @test isfile(validation.stats_path)
        @test isfile(validation.query_vs_uniprot_path)
        stats_row = only(collect(Iduna.ResultsValidation.CSV.File(validation.stats_path)))
        @test stats_row.seed_path == "tiny.sto"
        @test stats_row.expanded_path == "tiny.sto"
        @test stats_row.query_vs_uniprot_path ==
              joinpath("validation", "pid_10.00", "query_vs_uniprot_alignment.txt")

        indel_uniprot = joinpath(tmp, "uniprot_indel.fasta")
        write(indel_uniprot, ">P20963\nACD\n")
        indel_target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001",
            uniprot_sequence_path = indel_uniprot
        )
        indel_validation = Iduna.ResultsValidation.validate_results(
            indel_target, seed, expansion, joinpath(tmp, "indel"))
        @test indel_validation.status === :warn
        @test occursin("indels relative to the UniProt sequence",
            only(indel_validation.warnings))
        indel_row = (;
            aln_identical = false,
            aln_insertions = 1,
            aln_deletions = 0)
        @test occursin("Seed query has indels",
            only(Iduna.ResultsValidation._cached_alignment_warnings(indel_row, nothing)))
        identical_row = (;
            aln_identical = true,
            aln_insertions = 0,
            aln_deletions = 0)
        @test isempty(
            Iduna.ResultsValidation._cached_alignment_warnings(identical_row, nothing))
        @test isempty(Iduna.ResultsValidation._alignment_warnings(nothing,
            (; insertions = 0, deletions = 0, identical = true)))

        relative_seed = Iduna.SeedSelection(;
            pid = 10.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = "tiny.sto",
            summary_path = "seed_summary.csv",
            workdir = tmp)
        relative_expansion = Iduna.ExpansionResult(;
            run_dir = ".",
            seed_stockholm = "tiny.sto",
            hits_fasta = "tiny.fasta",
            full_stockholm = "tiny.sto",
            match_stockholm = "tiny.sto",
            a3m_path = "tiny.a3m",
            db_dir = ".",
            hmm_dir = ".",
            logs_dir = ".",
            workdir = tmp)
        relative_target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001",
            uniprot_sequence_path = "uniprot.fasta",
            workdir = tmp)
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_seed_msa(relative_seed)) == 2
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_expanded_msa(relative_expansion)) == 2
        relative_validation = Iduna.ResultsValidation.validate_results(
            relative_target, relative_seed, relative_expansion, tmp)
        @test relative_validation.status === :warn
        @test relative_validation.query_name == "seq1"

        seed_only_validation = Iduna.ResultsValidation.validate_results(
            target, seed, nothing, joinpath(tmp, "seed_only"))
        @test seed_only_validation.status === :warn
        @test seed_only_validation.query_name == "seq1"
        @test seed_only_validation.seed_nseq == 2
        @test seed_only_validation.expanded_nseq === nothing
        @test seed_only_validation.expanded_ncol === nothing
        @test seed_only_validation.expanded_clusters62 === nothing
        @test seed_only_validation.expanded_neff80 === nothing
        @test occursin("Seed query has substitutions", only(seed_only_validation.warnings))
        @test isfile(seed_only_validation.stats_path)
        @test isfile(seed_only_validation.query_vs_uniprot_path)

        transcript_target = Iduna.ResolvedTarget(;
            input_id = "ENST000001",
            input_kind = :ensembl_transcript,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001"
        )
        transcript_validation = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, expansion, joinpath(tmp, "transcript"))
        @test transcript_validation.status === :ok
        @test isfile(transcript_validation.stats_path)
        overwritten_transcript_validation = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, expansion, joinpath(tmp, "transcript");
            overwrite = true)
        @test overwritten_transcript_validation.status === :ok
        cached_transcript_validation = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, expansion, joinpath(tmp, "transcript"))
        @test cached_transcript_validation.status === :ok
        skipped_expansion = Iduna.ExpansionResult(;
            run_dir = expansion.run_dir,
            seed_stockholm = expansion.seed_stockholm,
            seed_fasta = expansion.seed_fasta,
            hits_fasta = expansion.hits_fasta,
            full_stockholm = expansion.full_stockholm,
            match_stockholm = expansion.match_stockholm,
            a3m_path = expansion.a3m_path,
            s_exon_blocks_tsv = expansion.s_exon_blocks_tsv,
            db_dir = expansion.db_dir,
            hmm_dir = expansion.hmm_dir,
            logs_dir = expansion.logs_dir,
            n_hits = expansion.n_hits,
            n_new_hits = expansion.n_new_hits,
            status = :skipped,
            workdir = expansion.workdir)
        cached_after_expansion_reuse = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, skipped_expansion, joinpath(tmp, "transcript"))
        @test cached_after_expansion_reuse.status === :ok
        validation_state = Iduna.Utils._read_stage_state(
            Iduna.ResultsValidation._validation_dir(joinpath(tmp, "transcript"), seed))
        @test validation_state["action"] == "reuse"

        missing_seed = Iduna.SeedSelection(;
            pid = 20.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = joinpath(tmp, "missing_seed.sto"),
            summary_path = joinpath(tmp, "missing_summary.csv"))
        failure_workdir = joinpath(tmp, "validation_failure")
        @test_throws ErrorException Iduna.ResultsValidation.validate_results(
            transcript_target, missing_seed, nothing, failure_workdir)
        failed_state = Iduna.Utils._read_stage_state(
            Iduna.ResultsValidation._validation_dir(failure_workdir, missing_seed))
        @test failed_state["status"] == "failed"
        @test failed_state["exception"]["type"] == "ErrorException"
    end
end
