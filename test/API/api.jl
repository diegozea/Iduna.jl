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

    best_seed = Iduna.SeedSelection(;
        pid = 10.0,
        median_identity = 1.0,
        mean_identity = 1.0,
        stockholm_path = "seed.sto",
        summary_path = "best_seed.csv")
    thoraxe = Iduna.ThorAxeMSAResult(;
        input_dir = "thoraxe_input",
        thoraxe_dir = "thoraxe",
        msa_dir = "thoraxe_msa",
        baseline_fasta = "msa_0.fasta",
        baseline_stockholm = "msa_0.sto",
        sequence_fasta = "msa_0_sequences.fasta",
        species_file = "species.txt",
        pid_summary = "best_seed.csv",
        best_seed,
        logs_dir = "logs/thoraxe",
        warnings = ["BioMart transcript_query failures recorded for species: mus_spretus."],
        status = :warn)
    target = Iduna.ResolvedTarget(;
        input_id = "Q13148",
        input_kind = :uniprot,
        uniprot_id = "Q13148",
        ensembl_gene_id = "ENSG00000120948.20",
        transcript_id = "ENST00000240185.8")
    expansion = Iduna.ExpansionResult(;
        run_dir = "expansion",
        seed_stockholm = "seed.sto",
        hits_fasta = "hits.fasta",
        full_stockholm = "full.sto",
        match_stockholm = "match.sto",
        a3m_path = "expanded.a3m",
        db_dir = "db",
        hmm_dir = "hmm",
        logs_dir = "logs/expansion")
    validation = Iduna.ValidationResult(; stats_path = "stats.csv")
    result = Iduna.IdunaResult(;
        input_id = "Q13148",
        workdir = "workdir",
        target,
        thoraxe_msa = thoraxe,
        expansion,
        validation,
        warnings = thoraxe.warnings,
        status = :warn)
    summary = Iduna.Utils.result_summary(result)
    @test summary.status == "warn"
    @test summary.thoraxe_msa.status == "warn"
    @test summary.thoraxe_msa.warnings == thoraxe.warnings

    @testset "centroids option forwarding" begin
        mktempdir() do tmp
            captured = Ref{Dict{Symbol, Any}}()
            result = Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "db",
                workdir = joinpath(tmp, "centroids_forwarding"),
                centroids = true,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                _expand_msa = (args...; kwargs...) -> begin
                    captured[] = Dict{Symbol, Any}(kwargs)
                    expansion
                end,
                _validate_results = (args...; kwargs...) -> validation)

            @test result.expansion === expansion
            @test captured[][:centroids] === true
            @test captured[][:mmseqs_db] == "db"
        end
    end

    @testset "failure result artifacts" begin
        _read_result(path) = Iduna.JSON3.read(read(joinpath(path, "result.json"), String))

        mktempdir() do tmp
            workdir = joinpath(tmp, "target_failure")
            target_failure = (args...; kwargs...) -> error("target boom")
            @test_throws ErrorException Iduna.iduna(;
                id = "P20963",
                mmseqs_db = "db",
                workdir,
                _resolve_target = target_failure)

            failed = _read_result(workdir)
            @test failed.status == "error"
            @test failed.failed_stage == "resolve_target"
            @test failed.target === nothing
            @test failed.exception.type == "ErrorException"
            @test failed.exception.message == "target boom"
            @test !isfile(joinpath(workdir, "target.json"))
        end

        mktempdir() do tmp
            workdir = joinpath(tmp, "thoraxe_failure")
            thoraxe_failure = (args...; kwargs...) -> error("thoraxe boom")
            @test_throws ErrorException Iduna.iduna(;
                id = "P20963",
                mmseqs_db = "db",
                workdir,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = thoraxe_failure)

            failed = _read_result(workdir)
            @test failed.status == "error"
            @test failed.failed_stage == "thoraxe_msa"
            @test failed.target.uniprot_id == "Q13148"
            @test failed.warnings == String[]
            @test failed.exception.type == "ErrorException"
            @test failed.exception.message == "thoraxe boom"
            @test isfile(joinpath(workdir, "target.json"))
        end

        mktempdir() do tmp
            workdir = joinpath(tmp, "timeout_failure")
            stdout_log = joinpath(workdir, "logs", "thoraxe", "stdout.log")
            stderr_log = joinpath(workdir, "logs", "thoraxe", "stderr.log")
            timeout_failure = (args...; kwargs...) -> throw(
                Iduna.ThorAxeMSA._CommandTimeoutError(
                    "thoraxe --example", 12.0, stdout_log, stderr_log))
            @test_throws Iduna.ThorAxeMSA._CommandTimeoutError Iduna.iduna(;
                id = "P20963",
                mmseqs_db = "db",
                workdir,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = timeout_failure)

            failed = _read_result(workdir)
            @test failed.status == "error"
            @test failed.failed_stage == "thoraxe_msa"
            @test failed.exception.type == "Iduna.ThorAxeMSA._CommandTimeoutError"
            @test failed.exception.command == "thoraxe --example"
            @test failed.exception.stdout_log == stdout_log
            @test failed.exception.stderr_log == stderr_log
            @test occursin("timed out after 12.0 seconds", failed.exception.message)
        end
    end
end
