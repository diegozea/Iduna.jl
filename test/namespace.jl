@testset "Namespace boundary" begin
    @test isdefined(Iduna, :iduna)
    @test isdefined(Iduna, :load_result)
    @test isdefined(Iduna, :load_seed_msa)
    @test isdefined(Iduna, :load_expanded_msa)
    @test isdefined(Iduna, :_pipeline_status)
    @test !isdefined(Iduna, :pipeline_status)
    @test !isdefined(Iduna, :normalize_primary_input)
    @test !isdefined(Iduna, :parse_app_args)

    @test isdefined(Iduna.Utils, :DEFAULT_PID_THRESHOLDS)
    @test isdefined(Iduna.Utils, :_PROTEIN_ALIGNMENT_SCORE_MODEL)
    @test !isdefined(Iduna.Utils, :PROTEIN_ALIGNMENT_SCORE_MODEL)
    @test !isdefined(Iduna.Utils, :wrap_sequence)

    @test isdefined(Iduna.IDMapping, :resolve_target)
    @test isdefined(Iduna.IDMapping, :_parse_transcript_lookup)
    @test !isdefined(Iduna.IDMapping, :parse_transcript_lookup)
    @test !isdefined(Iduna.IDMapping, :UNIPROT_BASE)

    @test isdefined(Iduna.ThorAxeMSA, :build_thoraxe_msa)
    @test !isdefined(Iduna.ThorAxeMSA, :REQUIRED_ENSEMBL_FILES)

    @test isdefined(Iduna.MSAExpansion, :expand_msa)
    @test isdefined(Iduna.MSAExpansion, :_collect_hits)
    @test !isdefined(Iduna.MSAExpansion, :collect_hits)

    @test isdefined(Iduna.ResultsValidation, :alignment_stats)
    @test isdefined(Iduna.ResultsValidation, :_resolve_query_name)
    @test !isdefined(Iduna.ResultsValidation, :resolve_query_name)
end
