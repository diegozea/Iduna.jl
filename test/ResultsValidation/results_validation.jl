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
        @test Iduna.ResultsValidation.resolve_query_name(msa, "missing_gene", "missing_tx") ==
              "seq1"
        stats = Iduna.ResultsValidation.alignment_stats(sto)
        @test stats.n_sequences == 2
        @test stats.n_columns == 4
    end
end
