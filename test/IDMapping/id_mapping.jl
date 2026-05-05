@testset "IDMapping" begin
    @test Iduna.IDMapping.sequences_match("ACDEFG", "ACDEFG")
    @test !Iduna.IDMapping.sequences_match("ACDEFG", "ACDFFG")

    lookup = Iduna.IDMapping._parse_transcript_lookup(
        Dict(
            "id" => "ENSMUST00000193812",
            "Parent" => "ENSMUSG00000102693",
            "species" => "mus_musculus"
        ),
        "ENSMUST00000193812.2"
    )
    @test lookup.transcript_id == "ENSMUST00000193812.2"
    @test lookup.ensembl_gene_id == "ENSMUSG00000102693"
    @test lookup.species == "mus_musculus"

    no_species = Iduna.IDMapping._parse_transcript_lookup(
        Dict("Parent" => "ENSMUSG00000102693"),
        "ENSMUST00000193812.2"
    )
    @test no_species.species === nothing
    @test_throws ErrorException Iduna.IDMapping._parse_transcript_lookup(
        Dict("species" => "mus_musculus"),
        "ENSMUST00000193812.2"
    )

    mktempdir() do tmp
        resolver = _ -> lookup
        target = Iduna.IDMapping.resolve_target(
            "ENSMUST00000193812.2";
            workdir = tmp,
            _transcript_metadata_resolver = resolver
        )
        @test target.ensembl_gene_id == "ENSMUSG00000102693"
        @test target.species == "mus_musculus"

        explicit_species = Iduna.IDMapping.resolve_target(
            "ENSMUST00000193812.2";
            workdir = tmp,
            species = "custom_species",
            _transcript_metadata_resolver = resolver
        )
        @test explicit_species.ensembl_gene_id == "ENSMUSG00000102693"
        @test explicit_species.species == "custom_species"

        explicit_gene = Iduna.IDMapping.resolve_target(
            "ENSMUST00000193812.2";
            workdir = tmp,
            ensembl_gene_id = "ENSMUSG00000000001",
            _transcript_metadata_resolver = resolver
        )
        @test explicit_gene.ensembl_gene_id == "ENSMUSG00000000001"
        @test explicit_gene.species == "mus_musculus"
    end
end
