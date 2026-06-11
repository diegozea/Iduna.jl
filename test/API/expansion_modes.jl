@testset "multiple seed expansion mode" begin
    mktempdir() do tmp
        seed1 = Iduna.SeedSelection(;
            pid = 10.0,
            epli = missing,
            stockholm_path = "seed10.sto",
            summary_path = "candidate_summary.csv")
        seed2 = Iduna.SeedSelection(;
            pid = 80.0,
            epli = missing,
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
        workdir = joinpath(tmp, "no_expansion")
        stale_expansion_dir = joinpath(workdir, "expansion", "old_gene",
            "old_transcript", "pid_80.00")
        stale_validation_dir = joinpath(workdir, "validation", "pid_80.00")
        Iduna.Utils._write_stage_state(stale_expansion_dir;
            stage = "msa_expansion",
            stage_key = "expansion:pid_80.00",
            status = :done,
            identity = (; stale = "expansion"),
            outputs = NamedTuple(),
            action = :run,
            workdir)
        Iduna.Utils._write_stage_state(stale_validation_dir;
            stage = "validation",
            stage_key = "validation:pid_80.00",
            status = :done,
            identity = (; stale = "validation"),
            outputs = NamedTuple(),
            action = :run,
            workdir)
        stale_expansion_hash = JSON.parse(
            read(joinpath(stale_expansion_dir, "stage_state.json"), String))["identity_hash"]
        stale_validation_hash = JSON.parse(
            read(joinpath(stale_validation_dir, "stage_state.json"), String))["identity_hash"]
        result = Iduna.iduna(;
            id = "Q13148",
            workdir,
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
                Iduna.Utils._write_stage_state(
                    Iduna.ResultsValidation._validation_dir(workdir, seed);
                    stage = "validation",
                    stage_key = "validation:$(Iduna.Utils.format_pid_dir(seed.pid))",
                    status = :done,
                    identity = (; current = seed.pid),
                    outputs = NamedTuple(),
                    action = :run,
                    workdir)
                validation
            end)

        @test expand_called[] === false
        @test validation_expansion[] === nothing
        @test isempty(result.expansions)
        @test isempty(Iduna.Utils.result_summary(result).expansions)
        written = JSON.parse(
            read(joinpath(result.workdir, "result.json"), String))
        @test isempty(written["expansions"])
        stage_keys = [stage.stage_key for stage in result.stages]
        @test "validation:pid_10.00" in stage_keys
        @test !("expansion:pid_80.00" in stage_keys)
        @test !("validation:pid_80.00" in stage_keys)
        written_stage_keys = [stage["stage_key"] for stage in written["stages"]]
        @test !("expansion:pid_80.00" in written_stage_keys)
        @test !("validation:pid_80.00" in written_stage_keys)
        result_state = JSON.parse(read(
            joinpath(
                workdir, ".iduna", "stages", "result", "stage_state.json"), String))
        @test !(stale_expansion_hash in result_state["identity"]["stage_hashes"])
        @test !(stale_validation_hash in result_state["identity"]["stage_hashes"])
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
