@testset "orphan ThorAxe candidates do not satisfy current stage identity" begin
    mktempdir() do tmp
        workdir = joinpath(tmp, "work")
        paths = Iduna.ThorAxeMSA._pid_sample_paths(workdir, 10.0, 0)
        mkpath(dirname(paths.fasta_path))
        write(paths.fasta_path, ">stale\nAA\n")
        write(paths.stockholm_path, "# STOCKHOLM 1.0\nstale AA\n//\n")
        summary_path = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir),
            "candidate_summary.csv")
        stage_identity = (; target = "current")
        stage_cache = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(
            workdir, summary_path, (;), stage_identity; overwrite = false)
        @test stage_cache.cache.status === :missing
        @test stage_cache.summary_matches === false
        prepared = Iduna.ThorAxeMSA._prepare_thoraxe_msa_stage!(
            workdir, summary_path, stage_identity, stage_cache; overwrite = false)
        @test prepared.local_artifacts_are_current === false
        @test prepared.force_pid_rerun === true
    end
end

@testset "stale ThorAxe manifest cleanup removes old MSA dir" begin
    mktempdir() do tmp
        workdir = joinpath(tmp, "work")
        msa_dir = Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir)
        mkpath(msa_dir)
        write(joinpath(msa_dir, "stale.txt"), "stale")
        summary_path = joinpath(msa_dir, "candidate_summary.csv")
        stage_identity = (; target = "current")
        stage_cache = (;
            cache = (; reusable = true, status = :done, warning = nothing),
            has_manifest = true,
            summary_matches = false)
        prepared = @test_logs (:warn, r"manifest matched") Iduna.ThorAxeMSA._prepare_thoraxe_msa_stage!(
            workdir, summary_path, stage_identity, stage_cache; overwrite = false)
        @test prepared.force_pid_rerun === true
        @test !isfile(joinpath(msa_dir, "stale.txt"))
    end
end

@testset "ThorAxe MSA stage failure records manifest" begin
    mktempdir() do tmp
        workdir = joinpath(tmp, "work")
        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1")
        input_dir = joinpath(tmp, "input")
        write_test_ensembl_bundle(input_dir)
        summary_path = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir),
            "candidate_summary.csv")
        stage_identity = (; target = target.ensembl_gene_id)
        filters = (;
            species_filter = (; warnings = String[]),
            biomart_filter = (; warnings = String[]))
        prepared = (; action = :run, force_pid_rerun = true)
        failing_runner = (args...; kwargs...) -> error("stage boom")
        @test_throws ErrorException Iduna.ThorAxeMSA._run_thoraxe_msa_stage_with_failure_state!(
            failing_runner, target, input_dir, workdir, summary_path, [10.0],
            nothing, (;), stage_identity, filters, prepared;
            pid_sample_count = 0,
            pid_sample_fraction = 0.8,
            sample_seed = UInt64(1),
            sampling_strategy = :independent)
        state = Iduna.Utils._read_stage_state(
            Iduna.ThorAxeMSA._thoraxe_msa_stage_dir(workdir))
        @test state["status"] == "failed"
        @test state["action"] == "run"
        @test state["exception"]["type"] == "ErrorException"
        @test state["exception"]["message"] == "stage boom"
    end
end
