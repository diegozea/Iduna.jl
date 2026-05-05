@testset "Utils" begin
    @test Iduna.Utils.id_kind("P20963") === :uniprot
    @test Iduna.Utils.id_kind("A0A0C4DGY3") === :uniprot
    @test Iduna.Utils.id_kind("A0A0C4DGY3-1") === :uniprot
    @test Iduna.Utils.id_kind("ENST00000362089.10") === :ensembl_transcript
    @test_throws ErrorException Iduna.Utils.id_kind("not-an-id")

    identical = Iduna.Utils.protein_alignment_stats("ACDE", "ACDE")
    @test identical.identical
    @test identical.mismatches == 0
    @test identical.insertions == 0
    @test identical.deletions == 0

    mismatch = Iduna.Utils.protein_alignment_stats("ACDF", "ACDE")
    @test !mismatch.identical
    @test mismatch.mismatches == 1

    indel = Iduna.Utils.protein_alignment_stats("ACDEF", "ACDE")
    @test !indel.identical
    @test indel.insertions + indel.deletions == 1

    mktempdir() do tmp
        sto = joinpath(tmp, "names.sto")
        write(sto, "# STOCKHOLM 1.0\n1|ENSG000001.13 AC\nENST00000362089 AC\n//\n")
        msa = Iduna.ResultsValidation.load_msa(sto)
        @test Iduna.Utils.resolve_sequence_name(msa, "ENSG000001.13") ==
              "1|ENSG000001.13"
        @test Iduna.Utils.resolve_sequence_name(msa, "ENSG000001") ==
              "1|ENSG000001.13"
        @test Iduna.Utils.resolve_sequence_name(msa, "ENST00000362089.10") ==
              "ENST00000362089"
        @test Iduna.Utils.resolve_sequence_name(msa, "missing") === nothing
        @test Iduna.Utils.resolve_sequence_name(msa, "missing"; fallback = true) ==
              "1|ENSG000001.13"

        workdir = Iduna.Utils.prepare_output_dir("P20963"; output_dir = joinpath(tmp, "P20963"))
        @test isdir(workdir)

        owned = joinpath(workdir, "owned")
        mkpath(owned)
        @test isdir(owned)
        Iduna.Utils.safe_rm(owned, workdir)
        @test !isdir(owned)
        @test_throws ErrorException Iduna.Utils.safe_rm(tmp, workdir)
    end
end
