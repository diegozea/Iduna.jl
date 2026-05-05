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
        @test Iduna.ThorAxeMSA.ncolumns(msa) == 4
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
end
