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

        @test_throws ErrorException Iduna.ThorAxeMSA._ensure_transcript_query(
            target, joinpath(tmp, "bad_work");
            cached_input_dir = joinpath(tmp, "missing"),
            overwrite = true)
    end

    @testset "orthology specieslist filter helpers" begin
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

        fallback = Iduna.ThorAxeMSA._resolve_effective_specieslist(
            target, "Canis lupus", "1:1";
            homology_species_fetcher = (target, orthology) -> error("temporary failure"))
        @test fallback.specieslist == "Canis lupus"
        @test length(fallback.warnings) == 1
        @test occursin("Ensembl specieslist filter failed", only(fallback.warnings))
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
