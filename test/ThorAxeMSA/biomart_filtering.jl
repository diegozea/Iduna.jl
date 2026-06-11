@testset "BioMart datasets filter" begin
    biomart_text = """

TableSet\thsapiens_gene_ensembl\tHuman genes (GRCh38.p14)\t1
TableSet\tmmusculus_gene_ensembl\tMouse genes (GRCm39)\t1
TableSet\tptroglodytes_gene_ensembl\tChimpanzee genes (Pan_tro_3.0)\t1
"""
    datasets = Iduna.ThorAxeMSA._parse_biomart_gene_datasets(biomart_text)
    @test "hsapiens_gene_ensembl" in datasets
    @test "mmusculus_gene_ensembl" in datasets
    @test "ptroglodytes_gene_ensembl" in datasets
    @test Iduna.ThorAxeMSA._biomart_gene_dataset_for_species("Homo sapiens") ==
          "hsapiens_gene_ensembl"
    @test Iduna.ThorAxeMSA._biomart_gene_dataset_for_species("mus_musculus") ==
          "mmusculus_gene_ensembl"

    target = Iduna.ResolvedTarget(;
        input_id = "ENST00000000001.1",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1",
        species = "Homo sapiens")
    loader = () -> (
        datasets = Set(["hsapiens_gene_ensembl", "mmusculus_gene_ensembl"]),
        warnings = String[])
    filtered = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
        target, "Homo sapiens,Mus musculus,Canis lupus";
        dataset_loader = loader)
    @test filtered.specieslist == "homo_sapiens,mus_musculus"
    @test length(filtered.warnings) == 1
    @test occursin("canis_lupus", only(filtered.warnings))

    missing_query = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
        target, "Homo sapiens,Mus musculus";
        dataset_loader = () -> (datasets = Set(["mmusculus_gene_ensembl"]),
            warnings = String[]))
    @test missing_query.specieslist == "homo_sapiens,mus_musculus"
    @test any(w -> occursin("query species homo_sapiens", w),
        missing_query.warnings)

    unchecked = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
        target, "homo_sapiens,canis_lupus_familiaris";
        dataset_loader = loader)
    @test unchecked.specieslist == "homo_sapiens,canis_lupus_familiaris"
    @test any(w -> occursin("possible species aliases", w), unchecked.warnings)

    fallback = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
        target, "Homo sapiens,Mus musculus";
        dataset_loader = () -> (datasets = nothing,
            warnings = ["BioMart datasets metadata refresh failed"]))
    @test fallback.specieslist == "Homo sapiens,Mus musculus"
    @test only(fallback.warnings) == "BioMart datasets metadata refresh failed"
end

@testset "BioMart datasets dated cache" begin
    response(status::Integer,
        body::AbstractString = "") = Iduna.Utils.HTTP.Response(status; body)

    @test isdir(Iduna.ThorAxeMSA._biomart_cache_dir())

    attempts = Ref(0)
    biomart_text = "TableSet\thsapiens_gene_ensembl\tHuman genes\t1\n"
    text = Iduna.ThorAxeMSA._fetch_biomart_datasets_text(;
        retries = 3,
        http_get = (url;
            kwargs...) -> begin
            attempts[] += 1
            attempts[] == 1 ? response(503) : response(200, biomart_text)
        end)
    @test text == biomart_text
    @test attempts[] == 2

    attempts[] = 0
    @test_throws ErrorException Iduna.ThorAxeMSA._fetch_biomart_datasets_text(;
        retries = 3,
        http_get = (url; kwargs...) -> begin
            attempts[] += 1
            response(404)
        end)
    @test attempts[] == 1

    mktempdir() do tmp
        biomart_text = "TableSet\thsapiens_gene_ensembl\tHuman genes\t1\n"
        fetches = Ref(0)
        fetcher = () -> begin
            fetches[] += 1
            biomart_text
        end
        loaded = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
            cache_dir = tmp,
            today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 6),
            fetcher)
        @test fetches[] == 1
        @test "hsapiens_gene_ensembl" in loaded.datasets
        @test isempty(loaded.warnings)

        same_day = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
            cache_dir = tmp,
            today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 6),
            fetcher = () -> error("should not refetch"))
        @test fetches[] == 1
        @test "hsapiens_gene_ensembl" in same_day.datasets

        stale = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
            cache_dir = tmp,
            today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 7),
            fetcher = () -> error("temporary failure"))
        @test "hsapiens_gene_ensembl" in stale.datasets
        @test length(stale.warnings) == 1
        @test occursin("stale cache from 2026-05-06", only(stale.warnings))
    end

    mktempdir() do tmp
        failed = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
            cache_dir = tmp,
            today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 6),
            fetcher = () -> error("temporary failure"))
        @test failed.datasets === nothing
        @test length(failed.warnings) == 1
        @test occursin("using the unfiltered specieslist", only(failed.warnings))
    end

    mktempdir() do tmp
        metadata_path = joinpath(tmp, Iduna.ThorAxeMSA._BIOMART_DATASETS_METADATA_FILE)
        @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) === nothing
        write(metadata_path, """{"download_date":"2026-05-06"}""")
        @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) == "2026-05-06"
        write(metadata_path, """{"download_time":"2026-05-06T12:00:00"}""")
        @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) === nothing
        write(metadata_path, "{not json")
        @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) === nothing
    end
end

@testset "BioMart transcript_query warnings" begin
    mktempdir() do tmp
        input_dir = joinpath(tmp, "thoraxe_input")
        ensembl_dir = joinpath(input_dir, "Ensembl")
        logs_dir = joinpath(tmp, "logs", "thoraxe")
        mkpath(ensembl_dir)
        mkpath(logs_dir)
        write(joinpath(ensembl_dir, "errors.csv"),
            "Species,GeneID\nmus_spretus,ENSMSPG00010016579\n")
        write(joinpath(logs_dir, "transcript_query_stderr.log"),
            """
            transcript_query.py:304: UserWarning: It can not found nomascus_leucogenys in biomart (tried: ['nleucogenys_gene_ensembl', 'nleucogenys_eg_gene']).
            Last response:
            Query ERROR: caught BioMart::Exception::Usage: Dataset nleucogenys_eg_gene NOT FOUND
              warnings.warn(...)
            transcript_query.py:728: UserWarning: Download failed for ENSPTRG00000006744 in pan_troglodytes!
              warnings.warn(...)
            transcript_query.py:728: UserWarning: Download failed for ENSMSPG00010016579 in mus_spretus!
              warnings.warn(...)
            """)
        warnings = Iduna.ThorAxeMSA._biomart_transcript_query_warnings(
            input_dir, logs_dir)
        @test length(warnings) == 1
        @test occursin("mus_spretus", only(warnings))
        @test occursin("nomascus_leucogenys", only(warnings))
        @test occursin("pan_troglodytes", only(warnings))
        @test occursin("errors.csv", only(warnings))
        @test occursin("transcript_query_stderr.log", only(warnings))

        malformed_errors = joinpath(input_dir, "Ensembl", "malformed_errors.csv")
        write(malformed_errors, "Other\nmus_spretus\n")
        @test isempty(Iduna.ThorAxeMSA._species_from_biomart_errors_file(
            malformed_errors))
    end
end
