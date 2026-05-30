import JSON

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
        summary_path = "candidate_summary.csv")
    thoraxe = Iduna.ThorAxeMSAResult(;
        input_dir = "thoraxe_input",
        thoraxe_dirs = ["thoraxe"],
        msa_dir = "thoraxe_msa",
        baseline_fastas = ["msa_0.fasta"],
        baseline_stockholms = ["msa_0.sto"],
        sequence_fastas = ["msa_0_sequences.fasta"],
        species_files = ["species.txt"],
        pid_summary = "candidate_summary.csv",
        seeds = [best_seed],
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
        expansions = [expansion],
        validations = [validation],
        warnings = thoraxe.warnings,
        status = :warn)
    summary = Iduna.Utils.result_summary(result)
    @test summary.status == "warn"
    @test summary.thoraxe_msa.status == "warn"
    @test summary.thoraxe_msa.warnings == thoraxe.warnings
    @test summary.expansions[1].match_stockholm == "match.sto"
    validation_warning = Iduna.ValidationResult(;
        stats_path = "warning_stats.csv",
        warnings = ["validation warning"])
    @test Iduna._partial_warnings(target, thoraxe, [validation_warning]) ==
          vcat(thoraxe.warnings, validation_warning.warnings)
    @test Iduna._select_seed_index(result; pid = nothing, index = 1, label = "seed") == 1
    @test Iduna._select_seed_index(result; pid = 10.0, index = nothing, label = "seed") == 1
    @test_throws ErrorException Iduna._select_seed_index(
        result; pid = nothing, index = 2, label = "seed")
    @test_throws ErrorException Iduna._select_seed_index(
        result; pid = 80.0, index = nothing, label = "seed")

    @testset "pipeline progress logging" begin
        mktempdir() do tmp
            logged = @test_logs (:info, r"Preparing Iduna output directory") (:info,
                r"Resolving target identifiers") (:info, r"Building ThorAxe MSA") (:info,
                r"Expanding MSA seed") (:info, r"Validating Iduna results") (:info,
                r"Writing Iduna result artifact") (:info,
                r"Iduna pipeline completed") match_mode=:any begin
                Iduna.iduna(;
                    id = "Q13148",
                    mmseqs_db = "db",
                    workdir = joinpath(tmp, "logged"),
                    _resolve_target = (args...; kwargs...) -> target,
                    _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                    _expand_msa = (args...; kwargs...) -> expansion,
                    _validate_results = (args...; kwargs...) -> validation)
            end
            @test logged.input_id == "Q13148"
            @test logged.status === :warn
        end
    end

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
        @test occursin(r"\n\s+expansions\s+=\s+\[ExpansionResult\(", result_text)
        @test occursin(r"\n\s+validations\s+=\s+\[ValidationResult\(", result_text)
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
                stockholm_path = joinpath(workdir, "thoraxe_msa", "candidates",
                    "pid_10.00", "candidate_msa_full.sto"),
                fasta_path = joinpath(workdir, "thoraxe_msa", "candidates", "pid_10.00",
                    "candidate_msa_full.fasta"),
                s_exon_blocks_tsv = joinpath(workdir, "thoraxe_msa", "candidates",
                    "pid_10.00", "candidate_msa_full_s_exon_blocks.tsv"),
                summary_path = joinpath(workdir, "thoraxe_msa", "candidate_summary.csv"))
            abs_thoraxe = Iduna.ThorAxeMSAResult(;
                input_dir = joinpath(workdir, "thoraxe_input"),
                thoraxe_dirs = [joinpath(workdir, "thoraxe")],
                msa_dir = joinpath(workdir, "thoraxe_msa"),
                baseline_fastas = [joinpath(workdir, "thoraxe_msa", "baseline.fasta")],
                baseline_stockholms = [joinpath(workdir, "thoraxe_msa", "baseline.sto")],
                sequence_fastas = [joinpath(workdir, "thoraxe_msa", "sequences.fasta")],
                species_files = [joinpath(workdir, "thoraxe_msa", "species.txt")],
                pid_summary = joinpath(workdir, "thoraxe_msa", "candidate_summary.csv"),
                seeds = [abs_seed],
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
                s_exon_blocks_tsv = joinpath(
                    expansion_dir, "expanded_msa", "$(transcript)_s_exon_blocks.tsv"),
                db_dir = joinpath(expansion_dir, "dbs"),
                hmm_dir = joinpath(expansion_dir, "hmm"),
                logs_dir = joinpath(expansion_dir, "logs"))
            abs_validation = Iduna.ValidationResult(;
                stats_path = joinpath(workdir, "validation", "pid_10.00", "stats.csv"),
                query_vs_uniprot_path = joinpath(
                    workdir, "validation", "pid_10.00", "query_vs_uniprot_alignment.txt"))
            mkpath(dirname(abs_seed.stockholm_path))
            write(abs_seed.stockholm_path, "# STOCKHOLM 1.0\nseed AC\n//\n")
            mkpath(dirname(abs_expansion.match_stockholm))
            write(abs_expansion.match_stockholm, "# STOCKHOLM 1.0\nseed AC\n//\n")

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
                  joinpath("thoraxe_msa", "candidate_summary.csv")
            @test result.thoraxe_msa.seeds[1].stockholm_path ==
                  joinpath("thoraxe_msa", "candidates", "pid_10.00",
                "candidate_msa_full.sto")
            @test result.thoraxe_msa.seeds[1].s_exon_blocks_tsv ==
                  joinpath("thoraxe_msa", "candidates", "pid_10.00",
                "candidate_msa_full_s_exon_blocks.tsv")
            @test result.thoraxe_msa.seeds[1].workdir == result.workdir
            @test result.expansions[1].match_stockholm == joinpath(
                "expansion", gene, transcript, "expanded_msa",
                "$(transcript)_matchonly.sto")
            @test result.expansions[1].s_exon_blocks_tsv == joinpath(
                "expansion", gene, transcript, "expanded_msa",
                "$(transcript)_s_exon_blocks.tsv")
            @test result.expansions[1].workdir == result.workdir
            @test result.validations[1].stats_path ==
                  joinpath("validation", "pid_10.00", "stats.csv")

            target_json = JSON.parse(read(joinpath(workdir, "target.json"), String))
            @test target_json["uniprot_sequence_path"] ==
                  joinpath("sequences", "uniprot", "Q13148.fasta")
            written = JSON.parse(read(joinpath(workdir, "result.json"), String))
            @test written["thoraxe_msa"]["pid_summary"] ==
                  joinpath("thoraxe_msa", "candidate_summary.csv")
            @test written["expansions"][1]["match_stockholm"] ==
                  result.expansions[1].match_stockholm
            @test written["expansions"][1]["s_exon_blocks_tsv"] ==
                  result.expansions[1].s_exon_blocks_tsv
            @test written["validations"][1]["stats_path"] ==
                  joinpath("validation", "pid_10.00", "stats.csv")
            @test Iduna.ResultsValidation.nsequences(Iduna.load_seed_msa(result)) == 1
            @test Iduna.ResultsValidation.nsequences(Iduna.load_expanded_msa(result)) == 1
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

            @test result.expansions[1].match_stockholm == expansion.match_stockholm
            @test captured[][:centroids] === true
            @test captured[][:mmseqs_db] == "db"
        end
    end

    @testset "PID sampling option forwarding" begin
        mktempdir() do tmp
            captured = Ref{Dict{Symbol, Any}}()
            result = Iduna.iduna(;
                id = "Q13148",
                workdir = joinpath(tmp, "pid_sampling_forwarding"),
                no_expansion = true,
                pid_sample_count = 12,
                pid_sample_fraction = 0.65,
                pid_sample_seed = 42,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (
                    args...; kwargs...) -> begin
                    captured[] = Dict{Symbol, Any}(kwargs)
                    thoraxe
                end,
                _validate_results = (args...; kwargs...) -> validation)

            @test isempty(result.expansions)
            @test captured[][:pid_sample_count] == 12
            @test captured[][:pid_sample_fraction] == 0.65
            @test captured[][:pid_sample_seed] == 42
        end
    end

    @testset "multiple seed expansion mode" begin
        mktempdir() do tmp
            seed1 = Iduna.SeedSelection(;
                pid = 10.0,
                median_identity = missing,
                mean_identity = missing,
                stockholm_path = "seed10.sto",
                summary_path = "candidate_summary.csv")
            seed2 = Iduna.SeedSelection(;
                pid = 80.0,
                median_identity = missing,
                mean_identity = missing,
                stockholm_path = "seed80.sto",
                summary_path = "candidate_summary.csv")
            multi_thoraxe = Iduna.ThorAxeMSAResult(;
                input_dir = "thoraxe_input",
                thoraxe_dirs = ["thoraxe10", "thoraxe80"],
                msa_dir = "thoraxe_msa",
                baseline_fastas = ["seed10.fasta", "seed80.fasta"],
                baseline_stockholms = ["seed10.sto", "seed80.sto"],
                sequence_fastas = ["seed10_sequences.fasta", "seed80_sequences.fasta"],
                species_files = ["seed10_species.txt", "seed80_species.txt"],
                pid_summary = "candidate_summary.csv",
                seeds = [seed1, seed2],
                logs_dir = "logs/thoraxe",
                pid_sample_count = 0)
            expanded_pids = Float64[]
            validated_pids = Float64[]
            result = Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "db",
                workdir = joinpath(tmp, "multi"),
                pid_sample_count = 0,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (args...; kwargs...) -> multi_thoraxe,
                _expand_msa = (target,
                    seed,
                    workdir;
                    kwargs...) -> begin
                    push!(expanded_pids, seed.pid)
                    Iduna.ExpansionResult(;
                        run_dir = "expansion/$(Iduna.Utils.format_pid_dir(seed.pid))",
                        seed_stockholm = seed.stockholm_path,
                        hits_fasta = "hits$(Int(seed.pid)).fasta",
                        full_stockholm = "full$(Int(seed.pid)).sto",
                        match_stockholm = "match$(Int(seed.pid)).sto",
                        a3m_path = "expanded$(Int(seed.pid)).a3m",
                        db_dir = "db",
                        hmm_dir = "hmm",
                        logs_dir = "logs")
                end,
                _validate_results = (target,
                    seed,
                    expansion_arg,
                    workdir) -> begin
                    push!(validated_pids, seed.pid)
                    @test expansion_arg !== nothing
                    Iduna.ValidationResult(;
                        stats_path = "validation/$(Iduna.Utils.format_pid_dir(seed.pid))/stats.csv")
                end)
            @test [seed.pid for seed in result.thoraxe_msa.seeds] == [10.0, 80.0]
            @test expanded_pids == [10.0, 80.0]
            @test validated_pids == [10.0, 80.0]
            @test length(result.expansions) == 2
            @test length(result.validations) == 2
            @test_throws ErrorException Iduna.load_seed_msa(result)
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
            @test isempty(result.expansions)
            @test isempty(Iduna.Utils.result_summary(result).expansions)
            written = JSON.parse(
                read(joinpath(result.workdir, "result.json"), String))
            @test isempty(written["expansions"])
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
            @test isempty(result.expansions)
        end
    end

    @testset "failure result artifacts" begin
        _read_result(path) = JSON.parse(read(joinpath(path, "result.json"), String))

        mktempdir() do tmp
            workdir = joinpath(tmp, "target_failure")
            target_failure = (args...; kwargs...) -> error("target boom")
            @test_throws ErrorException Iduna.iduna(;
                id = "P20963",
                mmseqs_db = "db",
                workdir,
                _resolve_target = target_failure)

            failed = _read_result(workdir)
            @test failed["status"] == "error"
            @test failed["failed_stage"] == "resolve_target"
            @test failed["target"] === nothing
            @test failed["exception"]["type"] == "ErrorException"
            @test failed["exception"]["message"] == "target boom"
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
            @test failed["status"] == "error"
            @test failed["failed_stage"] == "thoraxe_msa"
            @test failed["target"]["uniprot_id"] == "Q13148"
            @test failed["warnings"] == String[]
            @test failed["exception"]["type"] == "ErrorException"
            @test failed["exception"]["message"] == "thoraxe boom"
            @test isfile(joinpath(workdir, "target.json"))
        end
    end
end
