@testset "target and validation stage helpers" begin
    mktempdir() do tmp
        sto = joinpath(tmp, "seed.sto")
        write(sto, "# STOCKHOLM 1.0\nseq1 ACDE\n//\n")
        seed = Iduna.SeedSelection(;
            pid = 10.0,
            epli = 1.0,
            stockholm_path = sto,
            summary_path = joinpath(tmp, "summary.csv"))
        validation_target = Iduna.ResolvedTarget(;
            input_id = "seq1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "seq1",
            transcript_id = "seq1")
        called = Iduna._call_validate_results(
            Iduna.ResultsValidation.validate_results,
            validation_target, seed, nothing, tmp; overwrite = false)
        @test called.status === :ok

        unreadable_workdir = joinpath(tmp, "unreadable")
        mkpath(unreadable_workdir)
        write(joinpath(unreadable_workdir, "target.json"), "{bad json")
        @test Iduna._read_cached_target(unreadable_workdir) === nothing

        stale_workdir = joinpath(tmp, "stale")
        stale_target = Iduna.ResolvedTarget(;
            input_id = "Q13148",
            input_kind = :uniprot,
            uniprot_id = "Q13148",
            ensembl_gene_id = "ENSG00000120948.20",
            transcript_id = "ENST00000240185.8",
            uniprot_sequence_path = "missing.fasta")
        Iduna.Utils.write_json(joinpath(stale_workdir, "target.json"),
            Iduna._target_summary(stale_target, stale_workdir))
        @test Iduna._read_cached_target(stale_workdir) === nothing

        cleanup_workdir = joinpath(tmp, "cleanup")
        mkpath(joinpath(cleanup_workdir, "sequences"))
        write(joinpath(cleanup_workdir, "target.json"), "{}")
        Iduna._clear_target_outputs!(cleanup_workdir)
        @test !isfile(joinpath(cleanup_workdir, "target.json"))
        @test !isdir(joinpath(cleanup_workdir, "sequences"))

        blocked_workdir = joinpath(tmp, "blocked")
        write(blocked_workdir, "not a directory")
        @test_logs (:warn, r"Could not write failure result artifact") Iduna._write_failure_result(
            "Q13148", blocked_workdir, "result", ErrorException("boom"))
    end
end

@testset "pipeline logging uses caller logger" begin
    mktempdir() do tmp
        custom = Logging.ConsoleLogger(IOBuffer(), Logging.Warn)
        seen_logger = Ref{Any}()
        logged = Logging.with_logger(custom) do
            Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "db",
                workdir = joinpath(tmp, "logged"),
                _resolve_target = (args...;
                    kwargs...) -> begin
                    seen_logger[] = Logging.current_logger()
                    target
                end,
                _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                _expand_msa = (args...; kwargs...) -> expansion,
                _validate_results = (args...; kwargs...) -> validation)
        end
        @test logged.input_id == "Q13148"
        @test logged.status === :warn
        @test seen_logger[] === custom

        captured_workdir = joinpath(tmp, "captured")
        captured_logs,
        captured_result = Test.collect_test_logs() do
            Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "db",
                workdir = captured_workdir,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                _expand_msa = (args...; kwargs...) -> expansion,
                _validate_results = (args...; kwargs...) -> validation)
        end
        @test captured_result.status === :warn
        prepared_log = only([log
                             for log in captured_logs
                             if log.message == "Prepared Iduna output directory."])
        prepared_kwargs = Dict(prepared_log.kwargs)
        @test prepared_kwargs[:input_id] == "Q13148"
        @test prepared_kwargs[:workdir] == captured_workdir
        @test !haskey(prepared_kwargs, :overwrite)
        for message in ("Resolving target identifiers.",
            "Writing target metadata.",
            "Starting ThorAxe pipeline.",
            "Iduna pipeline completed.")
            log = only([log for log in captured_logs if log.message == message])
            kwargs = Dict(log.kwargs)
            @test !haskey(kwargs, :workdir)
            @test !haskey(kwargs, :input_id)
        end
        expanding_log = only([log
                              for log in captured_logs
                              if log.message == "Expanding MSA seed."])
        @test !haskey(Dict(expanding_log.kwargs), :pid)
        artifact_log = only([log
                             for log in captured_logs
                             if log.message == "Writing Iduna result artifact."])
        artifact_kwargs = Dict(artifact_log.kwargs)
        @test haskey(artifact_kwargs, :result_path)
        @test !haskey(artifact_kwargs, :status)

        overwrite_logs,
        _ = Test.collect_test_logs() do
            Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "db",
                workdir = joinpath(tmp, "captured_overwrite"),
                overwrite = true,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                _expand_msa = (args...; kwargs...) -> expansion,
                _validate_results = (args...; kwargs...) -> validation)
        end
        overwrite_prepared_log = only([log
                                       for log in overwrite_logs
                                       if log.message ==
                                          "Prepared Iduna output directory."])
        @test Dict(overwrite_prepared_log.kwargs)[:overwrite] === true
    end
end

@testset "result stage identity is stable across reruns" begin
    mktempdir() do tmp
        workdir = joinpath(tmp, "stable_result")
        for _ in 1:2
            Iduna.iduna(;
                id = "Q13148",
                mmseqs_db = "db",
                workdir,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
                _expand_msa = (args...; kwargs...) -> expansion,
                _validate_results = (args...; kwargs...) -> validation)
        end
        state_path = joinpath(workdir, ".iduna", "stages", "result",
            "stage_state.json")
        rerun_hash = JSON.parse(read(state_path, String))["identity_hash"]
        Iduna.iduna(;
            id = "Q13148",
            mmseqs_db = "db",
            workdir,
            _resolve_target = (args...; kwargs...) -> target,
            _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
            _expand_msa = (args...; kwargs...) -> expansion,
            _validate_results = (args...; kwargs...) -> validation)
        @test JSON.parse(read(state_path, String))["identity_hash"] == rerun_hash
        state = JSON.parse(read(state_path, String))
        @test state["output_digests"]["result.json"] ==
              Iduna.Utils._file_sha256(joinpath(workdir, "result.json"))
        @test Iduna.Utils._classify_stage_state(dirname(state_path), state["identity"],
            (; result_json = joinpath(workdir, "result.json"))).reusable
        saved = JSON.parse(read(joinpath(workdir, "result.json"), String))
        @test only(filter(stage -> stage["stage_key"] == "result", saved["stages"]))["status"] ==
              "done"
    end
end

@testset "target sequence files participate in cache integrity" begin
    mktempdir() do root
        calls = Ref(0)
        resolver = (primary; workdir, kwargs...) -> begin
            calls[] += 1
            sequence = joinpath(workdir, "sequences", "reference.fasta")
            mkpath(dirname(sequence))
            write(sequence, ">reference\nACDE\n")
            Iduna.ResolvedTarget(; input_id = primary, input_kind = :uniprot,
                ensembl_gene_id = "GENE", transcript_id = "TX",
                uniprot_sequence_path = sequence, workdir)
        end
        kwargs = (; supplied_uniprot = nothing, ensembl_gene_id = nothing,
            ensembl_protein_id = nothing, disambiguating_transcript = nothing,
            species = nothing, overwrite = false, resolver)
        first_target = Iduna._resolve_target_stage("P20963", root; kwargs...)
        Iduna._resolve_target_stage("P20963", root; kwargs...)
        @test calls[] == 1
        write(first_target.uniprot_sequence_path, ">reference\nACDF\n")
        rebuilt = Iduna._resolve_target_stage("P20963", root; kwargs...)
        @test calls[] == 2
        @test read(rebuilt.uniprot_sequence_path, String) == ">reference\nACDE\n"
        state_path = Iduna.Utils._stage_state_path(Iduna._target_stage_dir(root))
        state = JSON.parse(read(state_path, String))
        delete!(state, "output_digests")
        Iduna.Utils.write_json(state_path, state)
        Iduna._resolve_target_stage("P20963", root; kwargs...)
        @test calls[] == 3
        state = JSON.parse(read(state_path, String))
        delete!(state["required_outputs"], "uniprot_sequence")
        delete!(state["output_digests"], joinpath("sequences", "reference.fasta"))
        Iduna.Utils.write_json(state_path, state)
        Iduna._resolve_target_stage("P20963", root; kwargs...)
        @test calls[] == 4

        broken_resolver = (primary; workdir, kwargs...) -> Iduna.ResolvedTarget(;
            input_id = primary, input_kind = :uniprot,
            ensembl_gene_id = "GENE", transcript_id = "TX",
            uniprot_sequence_path = joinpath(workdir, "missing.fasta"))
        @test_throws ErrorException Iduna._resolve_target_stage("P20963", root;
            merge(kwargs, (; overwrite = true, resolver = broken_resolver))...)
    end
end

@testset "cached target stage preserves warnings" begin
    mktempdir() do tmp
        target_warning = "sequence and mapping warning"
        warned_target = Iduna.ResolvedTarget(;
            input_id = target.input_id,
            input_kind = target.input_kind,
            uniprot_id = target.uniprot_id,
            ensembl_gene_id = target.ensembl_gene_id,
            transcript_id = target.transcript_id,
            ensembl_protein_id = target.ensembl_protein_id,
            warnings = [target_warning])
        resolver_calls = Ref(0)
        workdir = joinpath(tmp, "target_warnings")
        run_kwargs = (;
            id = "Q13148",
            mmseqs_db = "db",
            workdir,
            _resolve_target = (args...; kwargs...) -> begin
                resolver_calls[] += 1
                warned_target
            end,
            _build_thoraxe_msa = (args...; kwargs...) -> thoraxe,
            _expand_msa = (args...; kwargs...) -> expansion,
            _validate_results = (args...; kwargs...) -> validation)

        Iduna.iduna(; run_kwargs...)
        rerun = Iduna.iduna(; run_kwargs...)

        @test resolver_calls[] == 1
        target_stage_path = joinpath(
            workdir, ".iduna", "stages", "target", "stage_state.json")
        target_state = JSON.parse(read(target_stage_path, String))
        @test target_state["action"] == "reuse"
        @test target_state["warnings"] == [target_warning]
        target_stages = filter(stage -> stage.stage_key == "target", rerun.stages)
        @test length(target_stages) == 1
        @test only(target_stages).warnings == [target_warning]
        written = JSON.parse(read(joinpath(workdir, "result.json"), String))
        written_target_stages = filter(
            stage -> stage["stage_key"] == "target", written["stages"])
        @test only(written_target_stages)["warnings"] == [target_warning]

        delete!(target_state["identity"], "target_contract")
        target_state["identity_hash"] = Iduna.Utils._identity_hash(target_state["identity"])
        Iduna.Utils.write_json(target_stage_path, target_state)
        Iduna.iduna(; run_kwargs...)
        @test resolver_calls[] == 2
        @test JSON.parse(read(target_stage_path, String))["action"] == "rebuild"
    end
end

@testset "target provenance controls reuse" begin
    mktempdir() do root
        calls = Ref(0)
        resolver = (primary; workdir, kwargs...) -> begin
            calls[] += 1
            Iduna.ResolvedTarget(; input_id = primary, input_kind = :uniprot,
                ensembl_gene_id = "GENE", transcript_id = "TX", workdir)
        end
        kwargs = (; supplied_uniprot = nothing, ensembl_gene_id = nothing,
            ensembl_protein_id = nothing, disambiguating_transcript = nothing,
            species = nothing, overwrite = false, resolver)
        Iduna._resolve_target_stage("P20963", root; kwargs...)
        state_path = Iduna.Utils._stage_state_path(Iduna._target_stage_dir(root))
        state = JSON.parse(read(state_path, String))
        @test state["identity"]["implementation"]["source_sha256"] !== nothing
        state["identity"]["implementation"]["source_sha256"] = "old-implementation"
        state["identity_hash"] = Iduna.Utils._identity_hash(state["identity"])
        Iduna.Utils.write_json(state_path, state)
        Iduna._resolve_target_stage("P20963", root; kwargs...)
        @test calls[] == 2
        Iduna._resolve_target_stage("P20963", root; kwargs...)
        @test calls[] == 2

        # A selected environment without the loaded dependencies cannot authorize reuse.
        previous_project = Base.active_project()
        try
            Iduna.Utils.Pkg.activate(joinpath(root, "empty-project"); io = devnull)
            for expected in (3, 4)
                Iduna._resolve_target_stage("P20963", root; kwargs...)
                @test calls[] == expected
                @test JSON.parse(read(state_path, String))["identity"]["implementation"] ===
                      nothing
            end
        finally
            Iduna.Utils.Pkg.activate(dirname(previous_project); io = devnull)
        end
        Iduna._resolve_target_stage("P20963", root; kwargs...)
        @test calls[] == 5
    end
end
