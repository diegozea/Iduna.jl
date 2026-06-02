@testset "IDMapping" begin
    response(body::AbstractString; status::Integer = 200) = Iduna.Utils.HTTP.Response(
        status; body)

    @test Iduna.IDMapping.sequences_match("ACDEFG", "ACDEFG")
    @test !Iduna.IDMapping.sequences_match("ACDEFG", "ACDFFG")
    @test Iduna.IDMapping._http_get("https://example.org/ok";
        http_get = (url; kwargs...) -> response("""{"ok":true}""")) !== nothing
    @test_logs (:warn, r"HTTP request failed") begin
        @test Iduna.IDMapping._http_get("https://example.org/retry";
            retries = 1,
            http_get = (url; kwargs...) -> response("busy"; status = 503)) === nothing
    end
    @test Iduna.IDMapping._http_get("https://example.org/missing";
        http_get = (url; kwargs...) -> response("missing"; status = 404)) === nothing

    uniprot_json = """
    {
      "organism": {"scientificName": "Homo sapiens"},
      "sequence": {"value": "acde"},
      "uniProtKBCrossReferences": [
        {"database": "PDB", "id": "1ABC"},
        {
          "database": "Ensembl",
          "id": "ENST000001.1",
          "isoformId": "P20963-1",
          "properties": [
            {"key": "GeneId", "value": "ENSG000001"},
            {"key": "ProteinId", "value": "ENSP000001.5"}
          ]
        },
        {
          "database": "Ensembl",
          "id": "ENST000002",
          "properties": [
            {"key": "ProteinId", "value": "ENSP000002"}
          ]
        },
        {
          "database": "Ensembl",
          "id": "ENST000001.1",
          "properties": [
            {"key": "GeneId", "value": "ENSG000001"}
          ]
        }
      ]
    }
    """
    entry = Iduna.IDMapping.fetch_uniprot_entry("P20963";
        _http_get_fn = url -> begin
            @test endswith(url, "/uniprotkb/P20963.json")
            response(uniprot_json)
        end)
    @test entry.id == "P20963"
    @test entry.species == "Homo sapiens"
    @test entry.gene_ids == ["ENSG000001"]
    @test entry.transcript_ids == ["ENST000001.1", "ENST000002"]
    @test entry.transcript_to_protein["ENST000001.1"] == "ENSP000001.5"
    @test entry.transcript_to_protein["ENST000002"] == "ENSP000002"
    @test entry.transcript_to_gene["ENST000002"] == "ENSG000001"
    @test entry.transcript_to_isoform["ENST000001.1"] == "P20963-1"
    @test entry.protein_sequence == "ACDE"
    @test_throws ErrorException Iduna.IDMapping.fetch_uniprot_entry(
        "P20963"; _http_get_fn = url -> nothing)

    @test Iduna.IDMapping._get_species(Dict()) === nothing
    @test Iduna.IDMapping._extract_uniprot_sequence(Dict()) === nothing
    @test Iduna.IDMapping._extract_uniprot_sequence(
        Dict("sequence" => Dict("value" => 42))) === nothing

    protein_seq = Iduna.IDMapping._fetch_ensembl_protein_sequence("ENSP000001.5";
        _http_get_fn = (url, headers;
            kwargs...) -> begin
            @test occursin("/sequence/id/ENSP000001?type=protein", url)
            @test headers == Iduna.IDMapping._FASTA_HEADERS
            response(">ENSP000001\nacde\n")
        end)
    @test protein_seq == "ACDE"
    @test Iduna.IDMapping._fetch_ensembl_protein_sequence("ENSP000001";
        _http_get_fn = (url, headers; kwargs...) -> nothing) === nothing

    @test Iduna.IDMapping._fetch_uniprot_fasta_sequence("P20963";
        _http_get_fn = (url, headers; retries,
            sleep_seconds) -> begin
            @test endswith(url, "/uniprotkb/P20963.fasta")
            @test headers == Iduna.IDMapping._FASTA_HEADERS
            @test retries == 5
            @test sleep_seconds == 2.0
            response(">P20963\nacde\n")
        end) == "ACDE"
    @test_throws ErrorException Iduna.IDMapping._fetch_uniprot_fasta_sequence(
        "P20963"; _http_get_fn = (url, headers; kwargs...) -> nothing)
    @test_throws ErrorException Iduna.IDMapping._fetch_uniprot_fasta_sequence(
        "P20963"; _http_get_fn = (url, headers; kwargs...) -> response(">empty\n"))

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
    fetched_lookup = Iduna.IDMapping._resolve_transcript_metadata(
        "ENSMUST00000193812.2";
        _http_get_fn = url -> begin
            @test occursin("/lookup/id/ENSMUST00000193812?expand=0", url)
            response("""{"Parent":"ENSMUSG00000102693","species":"mus_musculus"}""")
        end)
    @test fetched_lookup.ensembl_gene_id == "ENSMUSG00000102693"
    @test_throws ErrorException Iduna.IDMapping._resolve_transcript_metadata(
        "ENSMUST00000193812.2"; _http_get_fn = url -> nothing)

    mktempdir() do tmp
        candidate_entry = Iduna.IDMapping._UniProtEntry(;
            id = "P20963",
            species = "Homo sapiens",
            gene_ids = ["ENSG000001"],
            transcript_ids = ["ENST000001.1", "ENST000002", "ENST000003"],
            transcript_to_protein = Dict(
                "ENST000001.1" => "ENSP000001",
                "ENST000002" => "ENSP000002"
            ),
            transcript_to_gene = Dict(
                "ENST000001.1" => "ENSG000001",
                "ENST000002" => "ENSG000001"
            ),
            transcript_to_isoform = Dict("ENST000001.1" => "P20963-1"),
            protein_sequence = nothing
        )
        candidates,
        uniprot_path = Iduna.IDMapping._validate_candidates(
            candidate_entry,
            joinpath(tmp, "seqs");
            _uniprot_fasta_fetcher = id -> begin
                @test id == "P20963"
                "ACDE"
            end,
            _ensembl_protein_fetcher = protein -> protein == "ENSP000001" ? "ACDE" :
                                                  protein == "ENSP000002" ? "ACDF" :
                                                  nothing
        )
        @test length(candidates) == 1
        @test only(candidates).transcript_id == "ENST000001.1"
        @test only(candidates).ensembl_protein_id == "ENSP000001"
        @test only(candidates).isoform_id == "P20963-1"
        @test only(candidates).sequence_validated
        @test only(candidates).mapping_confirmed === true
        @test read(uniprot_path, String) == ">P20963\nACDE\n"
        @test isfile(joinpath(tmp, "seqs", "ensembl_proteins", "ENSP000001.fasta"))

        first_candidate = Iduna.IDMapping.EnsemblCandidate(;
            transcript_id = "ENST000001.1",
            ensembl_gene_id = "ENSG000001",
            ensembl_protein_id = "ENSP000001",
            species = "Homo sapiens",
            sequence_validated = true,
            mapping_confirmed = true
        )
        second_candidate = Iduna.IDMapping.EnsemblCandidate(;
            transcript_id = "ENST000002.4",
            ensembl_gene_id = "ENSG000001",
            ensembl_protein_id = "ENSP000002",
            species = "Homo sapiens",
            sequence_validated = true,
            mapping_confirmed = true
        )
        chosen,
        warnings = Iduna.IDMapping._choose_candidate(
            [first_candidate, second_candidate], nothing)
        @test chosen === first_candidate
        @test length(warnings) == 1
        @test occursin("Multiple validated", only(warnings))
        chosen,
        warnings = Iduna.IDMapping._choose_candidate(
            [first_candidate, second_candidate], "ENST000002")
        @test chosen === second_candidate
        @test isempty(warnings)
        @test_throws ErrorException Iduna.IDMapping._choose_candidate(
            Iduna.IDMapping.EnsemblCandidate[], nothing)
        @test_throws ErrorException Iduna.IDMapping._choose_candidate(
            [first_candidate], "ENST999999")

        resolved = Iduna.IDMapping.resolve_target(
            "P20963";
            workdir = tmp,
            ensembl_gene_id = "ENSG_OVERRIDE",
            ensembl_protein_id = "ENSP_OVERRIDE",
            species = "custom_species",
            _uniprot_entry_fetcher = id -> begin
                @test id == "P20963"
                candidate_entry
            end,
            _candidate_validator = (entry,
                sequence_dir) -> begin
                @test entry === candidate_entry
                @test endswith(sequence_dir, joinpath("sequences"))
                ([first_candidate], joinpath(sequence_dir, "uniprot", "P20963.fasta"))
            end
        )
        @test resolved.input_kind === :uniprot
        @test resolved.ensembl_gene_id == "ENSG_OVERRIDE"
        @test resolved.ensembl_protein_id == "ENSP_OVERRIDE"
        @test resolved.species == "custom_species"
        @test resolved.sequence_validated === true
        @test resolved.mapping_confirmed === true

        missing_gene_candidate = Iduna.IDMapping.EnsemblCandidate(;
            transcript_id = "ENST000004",
            ensembl_protein_id = "ENSP000004"
        )
        @test_throws ErrorException Iduna.IDMapping.resolve_target(
            "P20963";
            workdir = tmp,
            _uniprot_entry_fetcher = id -> candidate_entry,
            _candidate_validator = (entry,
                sequence_dir) -> (
                [missing_gene_candidate], joinpath(sequence_dir, "uniprot", "P20963.fasta"))
        )

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
