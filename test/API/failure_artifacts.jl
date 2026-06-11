@testset "current stage summaries distinguish expansion target" begin
    mktempdir() do tmp
        current_key = "expansion:$(target.ensembl_gene_id):$(target.transcript_id):pid_10.00"
        old_key = "expansion:OLDGENE:OLDTX:pid_10.00"
        legacy_key = "expansion:pid_10.00"
        for (key, dir) in (
            current_key => joinpath(tmp, "expansion", target.ensembl_gene_id,
                target.transcript_id, "pid_10.00"),
            old_key => joinpath(tmp, "expansion", "OLDGENE", "OLDTX", "pid_10.00"),
            legacy_key => joinpath(tmp, "expansion", "LEGACY", "LEGACYTX",
                "pid_10.00"))
            Iduna.Utils._write_stage_state(dir;
                stage = "msa_expansion",
                stage_key = key,
                status = :done,
                identity = (; key),
                outputs = NamedTuple(),
                action = :run,
                workdir = tmp)
        end
        keys = Iduna._current_stage_keys(target, thoraxe, [expansion], [validation])
        @test current_key in keys
        @test !(old_key in keys)
        @test !(legacy_key in keys)
        summaries = Iduna.Utils.collect_stage_summaries(tmp; stage_keys = keys)
        expansion_summaries = filter(
            summary -> summary.stage == "msa_expansion", summaries)
        @test length(expansion_summaries) == 1
        @test only(expansion_summaries).stage_key == current_key
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

    mktempdir() do tmp
        workdir = joinpath(tmp, "expansion_failure")
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
        expansion_failure = (target, seed, workdir;
            kwargs...) -> begin
            seed.pid == 80.0 && error("expansion boom")
            return expansion
        end
        @test_throws ErrorException Iduna.iduna(;
            id = "P20963",
            mmseqs_db = "db",
            workdir,
            pid_sample_count = 0,
            _resolve_target = (args...; kwargs...) -> target,
            _build_thoraxe_msa = (args...; kwargs...) -> multi_thoraxe,
            _expand_msa = expansion_failure)

        failed = _read_result(workdir)
        @test failed["status"] == "error"
        @test failed["failed_stage"] == "msa_expansion"
        @test length(failed["thoraxe_msa"]["seeds"]) == 2
        @test length(failed["expansions"]) == 2
        @test failed["expansions"][1]["match_stockholm"] == expansion.match_stockholm
        @test failed["expansions"][2] === nothing
    end
end
