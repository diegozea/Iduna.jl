@testset "orthology specieslist filter helpers" begin
    @test Iduna.ThorAxeMSA._normalize_species_name(nothing) === nothing
    @test Iduna.ThorAxeMSA._prepend_query_species(["mus_musculus"], nothing) ==
          ["mus_musculus"]
    @test Iduna.ThorAxeMSA._orthology_relationships("1:1") ==
          ["ortholog_one2one"]
    @test Iduna.ThorAxeMSA._orthology_relationships("1:n") ==
          ["ortholog_one2one", "ortholog_one2many"]
    @test Iduna.ThorAxeMSA._orthology_relationships("m:n") ==
          ["ortholog_one2one", "ortholog_one2many", "ortholog_many2many"]
    @test_throws ErrorException Iduna.ThorAxeMSA._orthology_relationships("all")

    homologies = [
        Dict(
            "type" => "ortholog_one2one",
            "target" => Dict("species" => "mus_musculus")),
        Dict(
            "type" => "ortholog_one2many",
            "target" => Dict("species" => "danio_rerio")),
        Dict(
            "type" => "ortholog_many2many",
            "target" => Dict("species" => "xenopus_tropicalis"))
    ]
    homology_data = Dict("data" => [Dict("homologies" => homologies)])
    @test Iduna.ThorAxeMSA._homology_species(homology_data, "1:1") ==
          ["mus_musculus"]
    @test Iduna.ThorAxeMSA._homology_species(homology_data, "1:n") ==
          ["mus_musculus", "danio_rerio"]
    @test Iduna.ThorAxeMSA._homology_species(homology_data, "m:n") ==
          ["mus_musculus", "danio_rerio", "xenopus_tropicalis"]

    unknown_species_target = Iduna.ResolvedTarget(;
        input_id = "ENST00000000001.1",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1")
    @test_throws ErrorException Iduna.ThorAxeMSA._fetch_ortholog_species(
        unknown_species_target, "1:1")

    species_target = Iduna.ResolvedTarget(;
        input_id = "ENST00000000001.1",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1",
        species = "Homo sapiens")
    fetched_species = Iduna.ThorAxeMSA._fetch_ortholog_species(
        species_target, "1:n";
        homology_data_fetcher = (
            species, gene_id) -> begin
            @test species == "homo_sapiens"
            @test gene_id == "ENSG00000000001.1"
            homology_data
        end)
    @test fetched_species == ["mus_musculus", "danio_rerio"]
end

@testset "species list parsing and filtering" begin
    mktempdir() do tmp
        species_file = joinpath(tmp, "species.txt")
        write(species_file, "Mus musculus\n\nDanio rerio\n")
        ases_specieslist = Iduna.ThorAxeMSA._ASES_DEFAULT_SPECIESLIST
        ases_resolved = Iduna.ThorAxeMSA._resolve_specieslist_preset("ases")
        @test ases_resolved.specieslist == ases_specieslist
        @test ases_resolved.mode === :ases
        @test Iduna.ThorAxeMSA._resolve_specieslist_preset(" ases ").specieslist ==
              ases_specieslist
        @test Iduna.ThorAxeMSA._resolve_specieslist_preset("all").specieslist ===
              nothing
        @test Iduna.ThorAxeMSA._resolve_specieslist_preset("").specieslist ===
              nothing
        @test Iduna.ThorAxeMSA._resolve_specieslist_preset("   ").specieslist ===
              nothing
        @test Iduna.ThorAxeMSA._resolve_specieslist_preset(
            "Homo sapiens, Mus musculus").specieslist ==
              "Homo sapiens, Mus musculus"
        @test Iduna.ThorAxeMSA._parse_specieslist(nothing) === nothing
        @test Iduna.ThorAxeMSA._parse_specieslist(
            "Homo sapiens,mus_musculus, Mus musculus ") ==
              ["homo_sapiens", "mus_musculus"]
        @test Iduna.ThorAxeMSA._parse_specieslist(ases_specieslist) ==
              Iduna.ThorAxeMSA._ASES_DEFAULT_SPECIES
        @test Iduna.ThorAxeMSA._parse_specieslist(species_file) ==
              ["mus_musculus", "danio_rerio"]
        @test Iduna.ThorAxeMSA._parse_specieslist("Canis lupus") ==
              ["canis_lupus"]

        ases_logs,
        _ = Test.collect_test_logs() do
            Iduna.ThorAxeMSA._log_specieslist_choice(ases_resolved)
        end
        ases_log = only(ases_logs)
        @test ases_log.message == "Using Ases default ThorAxe species list."
        ases_kwargs = Dict(ases_log.kwargs)
        @test ases_kwargs[:n_species] == 12
        @test ases_kwargs[:specieslist] == ases_specieslist

        all_logs,
        _ = Test.collect_test_logs() do
            Iduna.ThorAxeMSA._log_specieslist_choice(
                Iduna.ThorAxeMSA._resolve_specieslist_preset("all"))
        end
        @test only(all_logs).message ==
              "Using unrestricted ThorAxe species selection."
        empty_logs,
        _ = Test.collect_test_logs() do
            Iduna.ThorAxeMSA._log_specieslist_choice(
                Iduna.ThorAxeMSA._resolve_specieslist_preset(""))
        end
        @test only(empty_logs).message ==
              "Using unrestricted ThorAxe species selection."

        list_logs,
        _ = Test.collect_test_logs() do
            Iduna.ThorAxeMSA._log_specieslist_choice(
                Iduna.ThorAxeMSA._resolve_specieslist_preset(
                "Homo sapiens,Mus musculus"))
        end
        list_log = only(list_logs)
        @test list_log.message == "Using explicit ThorAxe species list."
        @test Dict(list_log.kwargs)[:n_species] == 2

        file_logs,
        _ = Test.collect_test_logs() do
            Iduna.ThorAxeMSA._log_specieslist_choice(
                Iduna.ThorAxeMSA._resolve_specieslist_preset(species_file))
        end
        file_log = only(file_logs)
        @test file_log.message == "Using ThorAxe species list file."
        file_kwargs = Dict(file_log.kwargs)
        @test file_kwargs[:specieslist_path] == species_file
        @test file_kwargs[:n_species] == 2

        single_logs,
        _ = Test.collect_test_logs() do
            Iduna.ThorAxeMSA._log_specieslist_choice(
                Iduna.ThorAxeMSA._resolve_specieslist_preset("Canis lupus"))
        end
        single_log = only(single_logs)
        @test single_log.message == "Using explicit ThorAxe species."
        @test Dict(single_log.kwargs)[:species] == "canis_lupus"
    end

    target = Iduna.ResolvedTarget(;
        input_id = "ENST00000000001.1",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1",
        species = "Homo sapiens")
    fetcher = (target, orthology) -> ["mus_musculus", "danio_rerio"]

    all_species = Iduna.ThorAxeMSA._resolve_effective_specieslist(
        target, nothing, "1:1"; homology_species_fetcher = fetcher)
    @test all_species.specieslist == "homo_sapiens,mus_musculus,danio_rerio"
    @test isempty(all_species.warnings)

    ases_default = Iduna.ThorAxeMSA._resolve_effective_specieslist(
        target, Iduna.ThorAxeMSA._ASES_DEFAULT_SPECIESLIST, "1:1";
        homology_species_fetcher = fetcher)
    @test ases_default.specieslist == "homo_sapiens,mus_musculus,danio_rerio"
    @test length(ases_default.warnings) == 1
    @test occursin("gorilla_gorilla", only(ases_default.warnings))

    filtered = Iduna.ThorAxeMSA._resolve_effective_specieslist(
        target, "Mus musculus,Canis lupus", "1:1";
        homology_species_fetcher = fetcher)
    @test filtered.specieslist == "homo_sapiens,mus_musculus"
    @test length(filtered.warnings) == 1
    @test occursin("canis_lupus", only(filtered.warnings))

    @test_throws ErrorException Iduna.ThorAxeMSA._resolve_effective_specieslist(
        target, "Canis lupus", "1:1"; homology_species_fetcher = fetcher)
    @test_throws ErrorException Iduna.ThorAxeMSA._resolve_effective_specieslist(
        target, nothing, "1:1";
        homology_species_fetcher = (target, orthology) -> String[])

    fallback = Iduna.ThorAxeMSA._resolve_effective_specieslist(
        target, "Canis lupus", "1:1";
        homology_species_fetcher = (target, orthology) -> error("temporary failure"))
    @test fallback.specieslist == "Canis lupus"
    @test length(fallback.warnings) == 1
    @test occursin("Ensembl specieslist filter failed", only(fallback.warnings))

    resolved_filters = Iduna.ThorAxeMSA._resolve_thoraxe_species_filters(
        target, "Mus musculus", "1:1", nothing, true, true;
        specieslist_resolver = (target,
            specieslist,
            orthology) -> (
            specieslist = "homo_sapiens,mus_musculus",
            warnings = ["species filter warning"]),
        biomart_resolver = (target,
            specieslist) -> (
            specieslist = specieslist,
            warnings = ["biomart filter warning"]))
    @test resolved_filters.effective_specieslist == "homo_sapiens,mus_musculus"
    @test resolved_filters.species_filter.warnings == ["species filter warning"]
    @test resolved_filters.biomart_filter.warnings == ["biomart filter warning"]
end
