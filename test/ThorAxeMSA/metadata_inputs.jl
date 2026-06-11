mktempdir() do tmp
    metadata_path = joinpath(tmp, "metadata.json")
    expected = (gene_id = "ENSG", specieslist = nothing)
    @test !Iduna.ThorAxeMSA._metadata_matches(metadata_path, expected)

    write(metadata_path, """{"gene_id":"ENSG","specieslist":null}""")
    @test Iduna.ThorAxeMSA._metadata_matches(metadata_path, expected)
    @test !Iduna.ThorAxeMSA._metadata_matches(
        metadata_path, (; gene_id = "OTHER", specieslist = nothing))
    @test !Iduna.ThorAxeMSA._metadata_matches(
        metadata_path, (; gene_id = "ENSG", transcript_id = "ENST"))

    write(metadata_path, "{not json")
    @test !Iduna.ThorAxeMSA._metadata_matches(metadata_path, expected)
    @test Iduna.ThorAxeMSA._retry_wait_seconds(6) == 30.0
end
