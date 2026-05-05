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
    end
end
