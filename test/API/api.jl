@testset "API" begin
    @test Iduna._pipeline_status(String[]) === :ok
    @test Iduna._pipeline_status(["target warning"]) === :warn
    @test Iduna._pipeline_status(["thoraxe warning"]) === :warn
    @test Iduna._pipeline_status(["validation warning"]) === :warn

    primary, transcript,
    supplied_uniprot = Iduna._normalize_primary_input(;
        id = "ENST00000362089.10",
        uniprot_id = "P20963",
        ensembl_transcript_id = nothing,
        transcript_id = nothing
    )
    @test primary == "ENST00000362089.10"
    @test transcript === nothing
    @test supplied_uniprot == "P20963"

    primary, transcript,
    supplied_uniprot = Iduna._normalize_primary_input(;
        id = nothing,
        uniprot_id = "P20963",
        ensembl_transcript_id = "ENST00000362089.10",
        transcript_id = nothing
    )
    @test primary == "ENST00000362089.10"
    @test transcript === nothing
    @test supplied_uniprot == "P20963"

    primary, transcript,
    supplied_uniprot = Iduna._normalize_primary_input(;
        id = nothing,
        uniprot_id = "P20963",
        ensembl_transcript_id = nothing,
        transcript_id = "ENST00000362089.10"
    )
    @test primary == "P20963"
    @test transcript == "ENST00000362089.10"
    @test supplied_uniprot == "P20963"

    primary, transcript,
    supplied_uniprot = Iduna._normalize_primary_input(;
        id = "P20963",
        uniprot_id = "P20963",
        ensembl_transcript_id = nothing,
        transcript_id = nothing
    )
    @test primary == "P20963"
    @test transcript === nothing
    @test supplied_uniprot == "P20963"

    primary, transcript,
    supplied_uniprot = Iduna._normalize_primary_input(;
        id = nothing,
        uniprot_id = nothing,
        ensembl_transcript_id = nothing,
        transcript_id = "ENST00000362089.10"
    )
    @test primary == "ENST00000362089.10"
    @test transcript === nothing
    @test supplied_uniprot === nothing

    @test_throws ErrorException Iduna._normalize_primary_input(;
        id = nothing,
        uniprot_id = nothing,
        ensembl_transcript_id = nothing,
        transcript_id = nothing
    )

    @test_throws ErrorException Iduna._normalize_primary_input(;
        id = "P20963",
        uniprot_id = "Q9UKL0",
        ensembl_transcript_id = nothing,
        transcript_id = nothing
    )

    @test_throws ErrorException Iduna._normalize_primary_input(;
        id = "ENST00000362089.10",
        uniprot_id = nothing,
        ensembl_transcript_id = "ENST00000392122.4",
        transcript_id = nothing
    )
end
