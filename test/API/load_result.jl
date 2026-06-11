@testset "load_result reconstructs current schema read-only" begin
    mktempdir() do tmp
        workdir = joinpath(tmp, "successful")
        fixture = _write_load_result_fixture(workdir)
        before = read.(
            (
                fixture.result_json,
                fixture.target_json,
                fixture.summary_path,
                fixture.stats_path,
                fixture.stage_state),
            String)

        loaded = Iduna.load_result(workdir)
        @test loaded isa Iduna.IdunaResult
        @test loaded.workdir == abspath(workdir)
        @test loaded.status === :ok
        @test loaded.target.ensembl_gene_id == fixture.target.ensembl_gene_id
        @test loaded.thoraxe_msa.seeds[1].pid == 10.0
        @test loaded.thoraxe_msa.pid_sample_seed == UInt64(7)
        @test loaded.thoraxe_msa.sampling_strategy === :independent
        @test loaded.expansions[1].n_hits == 0
        @test loaded.validations[1].query_name == "seed"
        @test Iduna.ResultsValidation.nsequences(Iduna.load_seed_msa(loaded)) == 1
        @test Iduna.ResultsValidation.nsequences(Iduna.load_expanded_msa(loaded)) == 1
        @test before == read.(
            (
                fixture.result_json,
                fixture.target_json,
                fixture.summary_path,
                fixture.stats_path,
                fixture.stage_state),
            String)
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "no_expansion")
        _write_load_result_fixture(workdir; expansion = false)
        loaded = Iduna.load_result(workdir)
        @test isempty(loaded.expansions)
        @test Iduna.ResultsValidation.nsequences(Iduna.load_seed_msa(loaded)) == 1
        @test_throws ErrorException Iduna.load_expanded_msa(loaded)
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "stage_summary")
        fixture = _write_load_result_fixture(workdir)
        data = JSON.parse(read(fixture.result_json, String))
        data["stages"] = Any[Dict(
            "stage" => "target",
            "stage_key" => "target",
            "status" => "done")]
        Iduna.Utils.write_json(fixture.result_json, data)

        loaded = Iduna.load_result(workdir)
        @test length(loaded.stages) == 1
        @test loaded.stages[1]["stage_key"] == "target"
    end

    mktempdir() do tmp
        original = joinpath(tmp, "original")
        moved = joinpath(tmp, "moved")
        _write_load_result_fixture(original)
        cp(original, moved)
        loaded = Iduna.load_result(moved)
        @test loaded.workdir == abspath(moved)
        @test loaded.thoraxe_msa.seeds[1].stockholm_path == joinpath(
            "thoraxe_msa", "candidates", "pid_10.00", "candidate_msa_full.sto")
        @test Iduna.ResultsValidation.nsequences(Iduna.load_seed_msa(loaded)) == 1
        @test Iduna.ResultsValidation.nsequences(Iduna.load_expanded_msa(loaded)) == 1
    end
end

@testset "load_result warns for partial objects" begin
    mktempdir() do tmp
        workdir = joinpath(tmp, "failed_after_target")
        fixture = _write_load_result_fixture(workdir)
        rm(fixture.target_json; force = true)
        rm(joinpath(workdir, "thoraxe_msa"); recursive = true, force = true)
        Iduna.Utils.write_json(joinpath(workdir, "result.json"),
            Iduna._failure_result_summary(
                "Q13148", workdir, "thoraxe_msa", ErrorException("boom");
                target = fixture.target))

        loaded = @test_logs (:warn, r"Loading partial IdunaResult") (:warn,
            r"Loading target metadata") (:warn,
            r"no ThorAxe MSA summary") match_mode=:any Iduna.load_result(workdir)
        @test loaded.status === :error
        @test loaded.target.transcript_id == fixture.target.transcript_id
        @test isempty(loaded.thoraxe_msa.seeds)
        @test isempty(loaded.expansions)
        @test isempty(loaded.validations)
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "failed_after_thoraxe")
        fixture = _write_load_result_fixture(workdir; validation = false,
            pids = [10.0, 80.0])
        summary = Iduna._failure_result_summary(
            "Q13148", workdir, "msa_expansion", ErrorException("boom");
            target = fixture.target,
            thoraxe = fixture.thoraxe,
            expansions = Union{Missing, Iduna.ExpansionResult}[fixture.expansions[1]])
        @test !(:workdir in propertynames(summary))
        @test summary.failed_stage == "msa_expansion"
        @test summary.thoraxe_msa.seeds[1].stockholm_path ==
              joinpath("thoraxe_msa", "candidates", "pid_10.00",
            "candidate_msa_full.sto")
        @test length(summary.expansions) == 2
        @test summary.expansions[1].match_stockholm ==
              relpath(fixture.expansions[1].match_stockholm, workdir)
        @test summary.expansions[2] === nothing
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "bad_target_json")
        fixture = _write_load_result_fixture(workdir; expansion = false,
            validation = false)
        write(fixture.target_json, "{bad json")
        loaded = @test_logs (:warn, r"Could not read target.json") (:warn,
            r"Loading target metadata") match_mode=:any Iduna.load_result(workdir)
        @test loaded.target.input_id == "Q13148"
        @test loaded.target.uniprot_sequence_path ==
              joinpath("sequences", "uniprot", "Q13148.fasta")
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "bad_target_summary")
        mkpath(workdir)
        Iduna.Utils.write_json(joinpath(workdir, "result.json"),
            Dict("input_id" => "Q13148", "status" => "error", "target" => Dict()))
        @test_logs (:warn, r"Loading partial IdunaResult") (:warn,
            r"Loading target metadata") (:warn,
            r"Could not reconstruct target metadata") match_mode=:any begin
            @test_throws ErrorException Iduna.load_result(workdir)
        end
    end

    mktempdir() do tmp
        @test_throws ErrorException Iduna.load_result(joinpath(tmp, "missing"))
        workdir = joinpath(tmp, "array_result")
        mkpath(workdir)
        write(joinpath(workdir, "result.json"), "[]")
        @test_throws ErrorException Iduna.load_result(workdir)
    end
end

@testset "load_result current schema partial summaries" begin
    mktempdir() do tmp
        workdir = joinpath(tmp, "unsupported_seed_schema")
        _write_load_result_fixture(workdir; expansion = false, validation = false)
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        delete!(data["thoraxe_msa"], "seeds")
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        @test_throws ErrorException Iduna.load_result(workdir)
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "bad_expansion")
        _write_load_result_fixture(workdir)
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        data["expansions"][1]["n_hits"] = "bad"
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = @test_logs (:warn, r"Could not reconstruct expansion") match_mode=:any Iduna.load_result(
            workdir)
        @test length(loaded.expansions) == 1
        @test ismissing(loaded.expansions[1])
        @test_throws ErrorException Iduna.load_expanded_msa(loaded)
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "default_expansion_paths")
        fixture = _write_load_result_fixture(workdir)
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        delete!(data["expansions"][1], "run_dir")
        delete!(data["expansions"][1], "match_stockholm")
        delete!(data["expansions"][1], "seed_fasta")
        delete!(data["expansions"][1], "s_exon_blocks_tsv")
        delete!(data["expansions"][1], "n_hits")
        delete!(data["expansions"][1], "n_new_hits")
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = Iduna.load_result(workdir)
        @test loaded.expansions[1].run_dir == joinpath("expansion",
            fixture.target.ensembl_gene_id,
            fixture.target.transcript_id,
            Iduna.Utils.format_pid_dir(fixture.seed.pid))
        @test loaded.expansions[1].seed_fasta == joinpath(
            loaded.expansions[1].run_dir, "seeds",
            "seed_pid$(Iduna.Utils.format_pid(fixture.seed.pid)).fasta")
        @test loaded.expansions[1].s_exon_blocks_tsv == joinpath(
            loaded.expansions[1].run_dir, "expanded_msa",
            "$(fixture.target.transcript_id)_s_exon_blocks.tsv")
        @test loaded.expansions[1].n_hits == 0
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "explicit_null_expansion_paths")
        _write_load_result_fixture(workdir)
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        data["expansions"][1]["seed_fasta"] = nothing
        data["expansions"][1]["s_exon_blocks_tsv"] = nothing
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = Iduna.load_result(workdir)
        @test loaded.expansions[1].seed_fasta === nothing
        @test loaded.expansions[1].s_exon_blocks_tsv === nothing
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "match_path_expansion_run_dir")
        fixture = _write_load_result_fixture(workdir)
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        delete!(data["expansions"][1], "run_dir")
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = Iduna.load_result(workdir)
        @test loaded.expansions[1].run_dir == relpath(
            dirname(dirname(fixture.expansions[1].match_stockholm)), workdir)
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "partial_multi_seed_expansion")
        _write_load_result_fixture(workdir; validation = false, pids = [10.0, 80.0])
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        data["expansions"][1]["n_hits"] = "bad"
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = @test_logs (:warn, r"Could not reconstruct expansion") match_mode=:any Iduna.load_result(
            workdir)
        @test length(loaded.expansions) == 2
        @test ismissing(loaded.expansions[1])
        @test loaded.expansions[2] isa Iduna.ExpansionResult
        @test_throws ErrorException Iduna.load_expanded_msa(loaded; index = 1)
        @test Iduna.ResultsValidation.nsequences(
            Iduna.load_expanded_msa(loaded; index = 2)) == 1
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "malformed_expansion_summaries")
        _write_load_result_fixture(workdir; validation = false,
            pids = [10.0, 80.0, 90.0])
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        data["expansions"][1] = nothing
        data["expansions"][2] = "partial"
        push!(data["expansions"], Dict("status" => "ok"))
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)

        loaded = @test_logs (:warn, r"Expansion summary is missing") (:warn,
            r"Skipping partial expansion summary") (:warn,
            r"Ignoring expansion summaries") match_mode=:any Iduna.load_result(workdir)
        @test length(loaded.expansions) == 3
        @test ismissing(loaded.expansions[1])
        @test ismissing(loaded.expansions[2])
        @test loaded.expansions[3] isa Iduna.ExpansionResult
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "failed_before_first_expansion")
        _write_load_result_fixture(workdir; expansion = false, validation = false,
            pids = [10.0, 80.0])
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        data["status"] = "error"
        data["failed_stage"] = "msa_expansion"
        data["exception"] = Dict("type" => "ErrorException", "message" => "boom")
        data["expansions"] = Any[]
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)

        loaded = @test_logs (:warn, r"Loading partial IdunaResult") (:warn,
            r"Expansion failed before producing summaries") match_mode=:any Iduna.load_result(
            workdir)
        @test length(loaded.expansions) == 2
        @test all(ismissing, loaded.expansions)
        err = try
            Iduna.load_expanded_msa(loaded; index = 1)
            nothing
        catch err
            err
        end
        @test err isa ErrorException
        @test occursin("No expanded MSA is available at index 1",
            sprint(showerror, err))
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "missing_validation_stats")
        fixture = _write_load_result_fixture(workdir; expansion = false)
        rm(fixture.stats_path; force = true)
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        merge!(data["validations"][1],
            Dict("query_name" => "fallback_seed",
                "seed_ncol" => 2,
                "seed_clusters62" => 1,
                "seed_neff80" => 1.0,
                "expanded_ncol" => 2,
                "expanded_clusters62" => 1,
                "expanded_neff80" => 1.0,
                "aln_identical" => true,
                "aln_mismatches" => 0,
                "aln_insertions" => 0,
                "aln_deletions" => 0,
                "warnings" => ["partial validation"]))
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = @test_logs (:warn, r"Validation stats are missing") match_mode=:any Iduna.load_result(
            workdir)
        @test loaded.validations[1].seed_ncol == 2
        @test loaded.validations[1].query_name == "fallback_seed"
        @test loaded.validations[1].warnings == ["partial validation"]
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "invalid_validation_stats")
        fixture = _write_load_result_fixture(workdir; expansion = false)
        write(fixture.stats_path, "bad\n1\n")
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        data["validations"][1]["seed_ncol"] = 2
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = @test_logs (:warn, r"Could not read validation stats") match_mode=:any Iduna.load_result(
            workdir)
        @test loaded.validations[1].seed_ncol == 2
    end

    mktempdir() do tmp
        workdir = joinpath(tmp, "malformed_validation_summaries")
        _write_load_result_fixture(workdir; expansion = false)
        data = JSON.parse(read(joinpath(workdir, "result.json"), String))
        data["validations"] = Any["partial", data["validations"][1]]
        Iduna.Utils.write_json(joinpath(workdir, "result.json"), data)
        loaded = @test_logs (:warn, r"Skipping partial validation summary") match_mode=:any Iduna.load_result(
            workdir)
        @test isempty(loaded.validations)
    end
end
