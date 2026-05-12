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
    @test summary.expansion.match_stockholm == "match.sto"

    @testset "result pretty printing" begin
        expansion_text = repr("text/plain", expansion)
        @test startswith(expansion_text, "ExpansionResult(\n")
        @test occursin(r"\n\s+run_dir\s+=\s+expansion", expansion_text)
        @test occursin(r"\n\s+n_hits\s+=\s+0", expansion_text)
        @test occursin(r"\n\s+status\s+=\s+:ok", expansion_text)

        validation_text = repr("text/plain", validation)
        @test startswith(validation_text, "ValidationResult(\n")
        @test occursin(r"\n\s+stats_path\s+=\s+stats\.csv", validation_text)
        @test occursin(r"\n\s+status\s+=\s+:ok", validation_text)

        result_text = repr("text/plain", result)
        @test startswith(result_text, "IdunaResult(\n")
        @test occursin(r"\n\s+target\s+=\s+ResolvedTarget\(", result_text)
        @test occursin(r"\n\s+expansion\s+=\s+ExpansionResult\(", result_text)
        @test occursin(r"\n\s+validation\s+=\s+ValidationResult\(", result_text)
    end

    @testset "returned result artifact paths are relative to workdir" begin
        mktempdir() do tmp
            workdir = joinpath(tmp, "relative_paths")
            gene = "ENSG00000120948.20"
            transcript = "ENST00000240185.8"
            abs_target = Iduna.ResolvedTarget(;
                input_id = "Q13148",
                input_kind = :uniprot,
                uniprot_id = "Q13148",
                ensembl_gene_id = gene,
                transcript_id = transcript,
                ensembl_protein_id = "ENSP00000240185",
                uniprot_sequence_path = joinpath(
                    workdir, "sequences", "uniprot", "Q13148.fasta"),
                ensembl_protein_sequence_path = joinpath(
                    workdir, "sequences", "ensembl_proteins", "ENSP00000240185.fasta"))
            abs_seed = Iduna.SeedSelection(;
                pid = 10.0,
                median_identity = 1.0,
                mean_identity = 1.0,
                stockholm_path = joinpath(workdir, "thoraxe_msa", "seeds", "pid10.sto"),
                fasta_path = joinpath(workdir, "thoraxe_msa", "seeds", "pid10.fasta"),
                summary_path = joinpath(workdir, "thoraxe_msa", "best_seed.csv"))
            abs_thoraxe = Iduna.ThorAxeMSAResult(;
                input_dir = joinpath(workdir, "thoraxe_input"),
                thoraxe_dir = joinpath(workdir, "thoraxe"),
                msa_dir = joinpath(workdir, "thoraxe_msa"),
                baseline_fasta = joinpath(workdir, "thoraxe_msa", "baseline.fasta"),
                baseline_stockholm = joinpath(workdir, "thoraxe_msa", "baseline.sto"),
                sequence_fasta = joinpath(workdir, "thoraxe_msa", "sequences.fasta"),
                species_file = joinpath(workdir, "thoraxe_msa", "species.txt"),
                pid_summary = joinpath(workdir, "thoraxe_msa", "best_seed.csv"),
                best_seed = abs_seed,
                logs_dir = joinpath(workdir, "logs", "thoraxe"))
            expansion_dir = joinpath(workdir, "expansion", gene, transcript)
            abs_expansion = Iduna.ExpansionResult(;
                run_dir = expansion_dir,
                seed_stockholm = joinpath(expansion_dir, "seeds", "seed.sto"),
                seed_fasta = joinpath(expansion_dir, "seeds", "seed.fasta"),
                hits_fasta = joinpath(
                    expansion_dir, "expanded_msa", "$(transcript)_hits_raw.fasta"),
                full_stockholm = joinpath(
                    expansion_dir, "expanded_msa", "$(transcript)_full.sto"),
                match_stockholm = joinpath(
                    expansion_dir, "expanded_msa", "$(transcript)_matchonly.sto"),
                a3m_path = joinpath(
                    expansion_dir, "expanded_msa", "$(transcript)_expanded.a3m"),
                db_dir = joinpath(expansion_dir, "dbs"),
                hmm_dir = joinpath(expansion_dir, "hmm"),
                logs_dir = joinpath(expansion_dir, "logs"))
            abs_validation = Iduna.ValidationResult(;
                stats_path = joinpath(workdir, "validation", "stats.csv"),
                query_vs_uniprot_path = joinpath(
                    workdir, "validation", "query_vs_uniprot_alignment.txt"))

            result = Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "db",
                workdir,
                _resolve_target = (args...; kwargs...) -> abs_target,
                _build_thoraxe_msa = (args...; kwargs...) -> abs_thoraxe,
                _expand_msa = (args...; kwargs...) -> abs_expansion,
                _validate_results = (args...; kwargs...) -> abs_validation)

            @test result.workdir == abspath(workdir)
            @test result.target.uniprot_sequence_path ==
                  joinpath("sequences", "uniprot", "Q13148.fasta")
            @test result.thoraxe_msa.pid_summary ==
                  joinpath("thoraxe_msa", "best_seed.csv")
            @test result.thoraxe_msa.best_seed.stockholm_path ==
                  joinpath("thoraxe_msa", "seeds", "pid10.sto")
            @test result.thoraxe_msa.best_seed.workdir == result.workdir
            @test result.expansion.match_stockholm == joinpath(
                "expansion", gene, transcript, "expanded_msa",
                "$(transcript)_matchonly.sto")
            @test result.expansion.workdir == result.workdir
            @test result.validation.stats_path == joinpath("validation", "stats.csv")

            target_json = Iduna.JSON3.read(read(joinpath(workdir, "target.json"), String))
            @test target_json.uniprot_sequence_path ==
                  joinpath("sequences", "uniprot", "Q13148.fasta")
            written = Iduna.JSON3.read(read(joinpath(workdir, "result.json"), String))
            @test written.thoraxe_msa.pid_summary ==
                  joinpath("thoraxe_msa", "best_seed.csv")
            @test written.expansion.match_stockholm == result.expansion.match_stockholm
            @test written.validation.stats_path == joinpath("validation", "stats.csv")
        end
    end

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
                _expand_msa = (
                    args...; kwargs...) -> begin
                    captured[] = Dict{Symbol, Any}(kwargs)
                    expansion
                end,
                _validate_results = (args...; kwargs...) -> validation)

            @test result.expansion.match_stockholm == expansion.match_stockholm
            @test captured[][:centroids] === true
            @test captured[][:mmseqs_db] == "db"
        end
    end

    @testset "no expansion option" begin
        @test_throws ErrorException Iduna.iduna(; id = "Q13148")
        @test_throws ErrorException Iduna.iduna(;
            id = "Q13148",
            no_expansion = true,
            centroids = true)

        mktempdir() do tmp
            expand_called = Ref(false)
            validation_expansion = Ref{Any}(:unset)
            result = Iduna.iduna(;
                id = "Q13148",
                workdir = joinpath(tmp, "no_expansion"),
                no_expansion = true,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                _expand_msa = (args...; kwargs...) -> begin
                    expand_called[] = true
                    expansion
                end,
                _validate_results = (target, seed, expansion_arg,
                    workdir) -> begin
                    validation_expansion[] = expansion_arg
                    validation
                end)

            @test expand_called[] === false
            @test validation_expansion[] === nothing
            @test result.expansion === nothing
            @test Iduna.Utils.result_summary(result).expansion === nothing
            written = Iduna.JSON3.read(
                read(joinpath(result.workdir, "result.json"), String))
            @test written.expansion === nothing
            @test_throws ErrorException Iduna.load_expanded_msa(result)
        end

        mktempdir() do tmp
            expand_called = Ref(false)
            result = Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "ignored_db",
                workdir = joinpath(tmp, "no_expansion_with_db"),
                no_expansion = true,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                _expand_msa = (args...; kwargs...) -> begin
                    expand_called[] = true
                    expansion
                end,
                _validate_results = (args...; kwargs...) -> validation)

            @test expand_called[] === false
            @test result.expansion === nothing
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
            timeout_failure = (args...;
                kwargs...) -> throw(
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
            @test failed.exception.stdout_log == joinpath("logs", "thoraxe", "stdout.log")
            @test failed.exception.stderr_log == joinpath("logs", "thoraxe", "stderr.log")
            @test occursin("timed out after 12.0 seconds", failed.exception.message)
        end
    end
end
