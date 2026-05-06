@testset "ThorAxeMSA" begin
    mktempdir() do tmp
        thoraxe = joinpath(tmp, "thoraxe")
        msa_dir = joinpath(thoraxe, "msa")
        mkpath(msa_dir)
        write(joinpath(thoraxe, "path_table.csv"),
            "TranscriptIDCluster,Path\nENST00000000001,start/1_0/2_0/stop\n")
        write(joinpath(thoraxe, "s_exon_table.csv"),
            "GeneID,Species\nENSG00000000001,homo_sapiens\nORTHO1,mus_musculus\n")
        write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
            ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")
        write(joinpath(msa_dir, "msa_s_exon_2_0.fasta"),
            ">ENSG00000000001\nCC\n>ORTHO1\nCD\n")

        msa,
        species = Iduna.ThorAxeMSA.assemble_transcript_msa(
            thoraxe, "ENSG00000000001.1", "ENST00000000001.1")
        @test Iduna.ThorAxeMSA.nsequences(msa) == 2
        name = first(Iduna.ThorAxeMSA.sequencenames(msa))
        @test length(Iduna.ThorAxeMSA.stringsequence(msa, name)) == 4
        @test species == ["homo_sapiens", "mus_musculus"]

        uniprot = joinpath(tmp, "uniprot.fasta")
        target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1",
            uniprot_sequence_path = uniprot)
        write(uniprot, ">P20963\nAACC\n")
        @test isempty(Iduna.ThorAxeMSA._validate_transcript_translation(target, msa))

        write(uniprot, ">P20963\nAAAC\n")
        warnings = Iduna.ThorAxeMSA._validate_transcript_translation(target, msa)
        @test length(warnings) == 1
        @test occursin("substitutions", only(warnings))

        write(uniprot, ">P20963\nAAC\n")
        @test_throws ErrorException Iduna.ThorAxeMSA._validate_transcript_translation(
            target, msa)

        missing_uniprot = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1",
            uniprot_sequence_path = joinpath(tmp, "missing.fasta"))
        missing_warnings = Iduna.ThorAxeMSA._validate_transcript_translation(
            missing_uniprot, msa)
        @test length(missing_warnings) == 1
        @test occursin("UniProt sequence file is missing", only(missing_warnings))

        no_uniprot = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1")
        @test isempty(Iduna.ThorAxeMSA._validate_transcript_translation(no_uniprot, msa))

        missing_reference = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "ENSG_MISSING",
            transcript_id = "ENST_MISSING",
            uniprot_sequence_path = uniprot)
        @test_throws ErrorException Iduna.ThorAxeMSA._validate_transcript_translation(
            missing_reference, msa)
    end

    mktempdir() do tmp
        summary = joinpath(tmp, "best_seed.csv")
        write(summary, """
        pid,median_identity,mean_identity,stockholm_path,fasta_path
        30.0,70.0,70.0,pid30.sto,pid30.fa
        10.0,70.0,80.0,pid10.sto,pid10.fa
        80.0,60.0,90.0,pid80.sto,pid80.fa
        """)
        seed = Iduna.ThorAxeMSA.select_best_seed(summary)
        @test seed.pid == 10.0
        @test seed.median_identity == 70.0
    end

    mktempdir() do tmp
        stdout_log = joinpath(tmp, "stdout.log")
        stderr_log = joinpath(tmp, "stderr.log")
        Iduna.ThorAxeMSA._run_logged_command(
            `$(Base.julia_cmd()) --startup-file=no -e "println(\"ok\")"`,
            stdout_log,
            stderr_log;
            timeout_seconds = 10)
        @test read(stdout_log, String) == "ok\n"
        @test isempty(read(stderr_log, String))

        no_timeout_stdout = joinpath(tmp, "no_timeout_stdout.log")
        no_timeout_stderr = joinpath(tmp, "no_timeout_stderr.log")
        Iduna.ThorAxeMSA._run_logged_command(
            `sh -c "printf no-timeout"`,
            no_timeout_stdout,
            no_timeout_stderr;
            timeout_seconds = nothing)
        @test read(no_timeout_stdout, String) == "no-timeout"
        @test isempty(read(no_timeout_stderr, String))

        @test_throws Iduna.ThorAxeMSA._CommandTimeoutError Iduna.ThorAxeMSA._run_logged_command(
            `$(Base.julia_cmd()) --startup-file=no -e "sleep(2)"`,
            joinpath(tmp, "slow_stdout.log"),
            joinpath(tmp, "slow_stderr.log");
            timeout_seconds = 0.1)
    end

    mktempdir() do tmp
        source = joinpath(tmp, "cached_input")
        ensembl = joinpath(source, "Ensembl")
        mkpath(ensembl)
        for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
            write(joinpath(ensembl, file), "x\n")
        end
        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1")
        workdir = joinpath(tmp, "work")
        copied = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
            cached_input_dir = source,
            overwrite = true)
        @test copied == joinpath(workdir, "thoraxe_input")
        @test isfile(joinpath(copied, "Ensembl", "sequences.fasta"))
        @test isdir(source)
        reused = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
            cached_input_dir = source,
            overwrite = false)
        @test reused == copied

        @test_throws ErrorException Iduna.ThorAxeMSA._ensure_transcript_query(
            target, joinpath(tmp, "bad_work");
            cached_input_dir = joinpath(tmp, "missing"),
            overwrite = true)
    end

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
    end

    @testset "species list parsing and filtering" begin
        mktempdir() do tmp
            species_file = joinpath(tmp, "species.txt")
            write(species_file, "Mus musculus\n\nDanio rerio\n")
            @test Iduna.ThorAxeMSA._parse_specieslist(nothing) === nothing
            @test Iduna.ThorAxeMSA._parse_specieslist(
                "Homo sapiens,mus_musculus, Mus musculus ") ==
                  ["homo_sapiens", "mus_musculus"]
            @test Iduna.ThorAxeMSA._parse_specieslist(species_file) ==
                  ["mus_musculus", "danio_rerio"]
            @test Iduna.ThorAxeMSA._parse_specieslist("Canis lupus") ==
                  ["canis_lupus"]
        end

        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1",
            species = "Homo sapiens")
        fetcher = (target, orthology) -> ["mus_musculus", "danio_rerio"]

        default = Iduna.ThorAxeMSA._resolve_effective_specieslist(
            target, nothing, "1:1"; homology_species_fetcher = fetcher)
        @test default.specieslist == "homo_sapiens,mus_musculus,danio_rerio"
        @test isempty(default.warnings)

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
    end

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
        response(status::Integer, body::AbstractString = "") =
            Iduna.ThorAxeMSA.HTTP.Response(status, Vector{UInt8}(body))

        attempts = Ref(0)
        biomart_text = "TableSet\thsapiens_gene_ensembl\tHuman genes\t1\n"
        text = Iduna.ThorAxeMSA._fetch_biomart_datasets_text(;
            retries = 3,
            http_get = (url; kwargs...) -> begin
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
        end
    end

    @testset "Ensembl homology download retries" begin
        response(status::Integer, body::AbstractString = "") =
            Iduna.ThorAxeMSA.HTTP.Response(status, Vector{UInt8}(body))

        attempts = Ref(0)
        data = Iduna.ThorAxeMSA._fetch_ensembl_homology_data(
            "homo_sapiens", "ENSG00000000001.1";
            retries = 3,
            sleep_seconds = 0,
            http_get = (url; headers, retry, status_exception) -> begin
                attempts[] += 1
                @test occursin("/homology/id/homo_sapiens/ENSG00000000001", url)
                @test headers == Iduna.ThorAxeMSA._ENSEMBL_JSON_HEADERS
                @test retry == false
                @test status_exception == false
                attempts[] == 1 ? response(500) : response(200, """{"data": []}""")
            end)
        @test get(data, :data, nothing) !== nothing
        @test attempts[] == 2

        attempts[] = 0
        @test_throws ErrorException Iduna.ThorAxeMSA._fetch_ensembl_homology_data(
            "homo_sapiens", "ENSG00000000001.1";
            retries = 3,
            sleep_seconds = 0,
            http_get = (url; kwargs...) -> begin
                attempts[] += 1
                response(400)
            end)
        @test attempts[] == 1
    end

    @testset "transcript_query receives orthology and effective species list" begin
        mktempdir() do tmp
            captured = Ref{Cmd}()
            runner = command -> (captured[] = command)
            Iduna.ThorAxeMSA._run_transcript_query_once(
                "ENSG00000000001", tmp, "homo_sapiens",
                "homo_sapiens,mus_musculus",
                joinpath(tmp, "stdout.log"), joinpath(tmp, "stderr.log");
                timeout_seconds = nothing,
                orthology = "1:n",
                runner = runner)
            parts = captured[].exec
            @test "--orthology" in parts
            @test parts[findfirst(==("--orthology"), parts) + 1] == "1:n"
            @test "--specieslist" in parts
            @test parts[findfirst(==("--specieslist"), parts) + 1] ==
                  "homo_sapiens,mus_musculus"
        end
    end
end
