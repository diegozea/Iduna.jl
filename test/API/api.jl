import JSON
import Logging

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
        epli = 1.0,
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
    validation = Iduna.ValidationResult(;
        stats_path = "stats.csv",
        seed_nseq = 100,
        seed_ncol = 1000,
        expanded_nseq = 250,
        expanded_ncol = 1000)
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
    @test !(:workdir in propertynames(summary))
    @test summary.status == "warn"
    @test summary.thoraxe_msa.status == "warn"
    @test summary.thoraxe_msa.warnings == thoraxe.warnings
    @test summary.thoraxe_msa.sampling_strategy == "independent"
    @test summary.thoraxe_msa.seeds[1].stockholm_path == "seed.sto"
    @test summary.expansions[1].match_stockholm == "match.sto"
    missing_expansion_result = Iduna.IdunaResult(;
        input_id = "Q13148",
        workdir = "workdir",
        target,
        thoraxe_msa = thoraxe,
        expansions = Union{Missing, Iduna.ExpansionResult}[missing],
        validations = Iduna.ValidationResult[],
        status = :error)
    @test Iduna.Utils.result_summary(missing_expansion_result).expansions[1] === nothing
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

    function _write_load_result_fixture(workdir;
            expansion::Bool = true,
            validation::Bool = true,
            pids = [10.0])
        gene = "ENSG00000120948.20"
        transcript = "ENST00000240185.8"
        fixture_target = Iduna.ResolvedTarget(;
            input_id = "Q13148",
            input_kind = :uniprot,
            uniprot_id = "Q13148",
            ensembl_gene_id = gene,
            transcript_id = transcript,
            ensembl_protein_id = "ENSP00000240185",
            uniprot_sequence_path = joinpath(
                workdir, "sequences", "uniprot", "Q13148.fasta"))
        mkpath(dirname(fixture_target.uniprot_sequence_path))
        write(fixture_target.uniprot_sequence_path, ">Q13148\nAC\n")
        Iduna.Utils.write_json(joinpath(workdir, "target.json"),
            Iduna._target_summary(fixture_target, workdir))

        summary_path = joinpath(workdir, "thoraxe_msa", "candidate_summary.csv")
        mkpath(dirname(summary_path))
        summary_rows = String[]
        fixture_seeds = Iduna.SeedSelection[]
        for pid in pids
            pid_dir = Iduna.Utils.format_pid_dir(pid)
            seed_sto = joinpath(workdir, "thoraxe_msa", "candidates", pid_dir,
                "candidate_msa_full.sto")
            seed_fasta = joinpath(workdir, "thoraxe_msa", "candidates", pid_dir,
                "candidate_msa_full.fasta")
            seed_blocks = joinpath(workdir, "thoraxe_msa", "candidates", pid_dir,
                "candidate_msa_full_s_exon_blocks.tsv")
            sequence_fasta = joinpath(workdir, "thoraxe_msa", "candidates", pid_dir,
                "sequences", "candidate_sequences_full.fasta")
            species_file = joinpath(workdir, "thoraxe_msa", "candidates", pid_dir,
                "species", "candidate_species_full.txt")
            mkpath(dirname(seed_sto))
            mkpath(dirname(sequence_fasta))
            mkpath(dirname(species_file))
            write(seed_sto, "# STOCKHOLM 1.0\nseed AC\n//\n")
            write(seed_fasta, ">seed\nAC\n")
            write(seed_blocks,
                "alignment\tpid\tcode\ts_exon_id\tstart_col\tend_col\tn_columns\n")
            write(sequence_fasta, ">seed\nAC\n")
            write(species_file, "homo_sapiens\n")
            push!(summary_rows,
                "$(pid),true,1.0,$(relpath(seed_sto, workdir)),$(relpath(seed_fasta, workdir)),$(relpath(sequence_fasta, workdir)),$(relpath(species_file, workdir))")
            push!(fixture_seeds,
                Iduna.SeedSelection(;
                    pid = Float64(pid),
                    epli = 1.0,
                    stockholm_path = seed_sto,
                    fasta_path = seed_fasta,
                    s_exon_blocks_tsv = seed_blocks,
                    summary_path))
        end
        write(summary_path,
            "pid,selected,epli,stockholm_path,fasta_path,sequence_fasta,species_file\n" *
            join(summary_rows, "\n") * "\n")
        fixture_thoraxe = Iduna.ThorAxeMSAResult(;
            input_dir = joinpath(workdir, "thoraxe_input"),
            thoraxe_dirs = [joinpath(workdir, "thoraxe_msa", "runs",
                                Iduna.Utils.format_pid_dir(seed.pid), "full", "thoraxe")
                            for seed in fixture_seeds],
            msa_dir = joinpath(workdir, "thoraxe_msa"),
            baseline_fastas = [seed.fasta_path for seed in fixture_seeds],
            baseline_stockholms = [seed.stockholm_path for seed in fixture_seeds],
            sequence_fastas = [joinpath(workdir, "thoraxe_msa", "candidates",
                                   Iduna.Utils.format_pid_dir(seed.pid), "sequences",
                                   "candidate_sequences_full.fasta")
                               for seed in fixture_seeds],
            species_files = [joinpath(workdir, "thoraxe_msa", "candidates",
                                 Iduna.Utils.format_pid_dir(seed.pid), "species",
                                 "candidate_species_full.txt")
                             for seed in fixture_seeds],
            pid_summary = summary_path,
            seeds = fixture_seeds,
            logs_dir = joinpath(workdir, "logs", "thoraxe"),
            pid_sample_count = 1,
            pid_sample_fraction = 1.0,
            pid_sample_seed = UInt64(7))

        fixture_expansions = Union{Missing, Iduna.ExpansionResult}[]
        if expansion
            for seed in fixture_seeds
                pid_dir = Iduna.Utils.format_pid_dir(seed.pid)
                run_dir = joinpath(workdir, "expansion", gene, transcript, pid_dir)
                expanded_dir = joinpath(run_dir, "expanded_msa")
                mkpath(expanded_dir)
                match_sto = joinpath(expanded_dir, "$(transcript)_matchonly.sto")
                full_sto = joinpath(expanded_dir, "$(transcript)_full.sto")
                hits_fasta = joinpath(expanded_dir, "$(transcript)_hits_raw.fasta")
                write(match_sto, "# STOCKHOLM 1.0\nseed AC\n//\n")
                write(full_sto, "# STOCKHOLM 1.0\nseed AC\n//\n")
                write(hits_fasta, "")
                push!(fixture_expansions,
                    Iduna.ExpansionResult(;
                        run_dir,
                        seed_stockholm = joinpath(
                            run_dir, "seeds", "seed_pid$(Iduna.Utils.format_pid(seed.pid)).sto"),
                        seed_fasta = joinpath(
                            run_dir, "seeds", "seed_pid$(Iduna.Utils.format_pid(seed.pid)).fasta"),
                        hits_fasta,
                        full_stockholm = full_sto,
                        match_stockholm = match_sto,
                        a3m_path = joinpath(expanded_dir, "$(transcript)_expanded.a3m"),
                        s_exon_blocks_tsv = joinpath(
                            expanded_dir, "$(transcript)_s_exon_blocks.tsv"),
                        db_dir = joinpath(run_dir, "dbs"),
                        hmm_dir = joinpath(run_dir, "hmm"),
                        logs_dir = joinpath(run_dir, "logs"),
                        n_hits = 0,
                        n_new_hits = 0))
            end
        end

        fixture_validations = Iduna.ValidationResult[]
        if validation
            for seed in fixture_seeds
                pid_dir = Iduna.Utils.format_pid_dir(seed.pid)
                stats_path = joinpath(workdir, "validation", pid_dir, "stats.csv")
                mkpath(dirname(stats_path))
                write(stats_path,
                    "query_name,seed_nseq,seed_ncol,seed_clusters62,seed_neff80,expanded_nseq,expanded_ncol,expanded_clusters62,expanded_neff80,aln_identical,aln_mismatches,aln_insertions,aln_deletions\nseed,1,2,1,1.0,1,2,1,1.0,true,0,0,0\n")
                push!(fixture_validations, Iduna.ValidationResult(; stats_path))
            end
        end

        fixture_result = Iduna.IdunaResult(;
            input_id = "Q13148",
            workdir,
            target = fixture_target,
            thoraxe_msa = fixture_thoraxe,
            expansions = fixture_expansions,
            validations = fixture_validations,
            status = :ok)
        Iduna.Utils.write_json(joinpath(workdir, "result.json"),
            Iduna.Utils.result_summary(fixture_result))
        stage_state = joinpath(workdir, ".iduna", "stages", "result",
            "stage_state.json")
        mkpath(dirname(stage_state))
        write(stage_state, "{\"stage\":\"result\"}")
        stats_path = validation ? fixture_validations[1].stats_path : nothing
        return (;
            target = fixture_target,
            thoraxe = fixture_thoraxe,
            seed = first(fixture_seeds),
            seeds = fixture_seeds,
            expansions = fixture_expansions,
            summary_path,
            stats_path,
            result_json = joinpath(workdir, "result.json"),
            target_json = joinpath(workdir, "target.json"),
            stage_state)
    end

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
        @test startswith(result_text, "IdunaResult Q13148 [:warn]\n\n")
        @test occursin(r"FIELD\s+SUMMARY", result_text)
        for field in (
            "input_id", "workdir", "target", "thoraxe_msa", "expansions", "validations",
            "stages", "warnings", "status")
            @test occursin(field, result_text)
            @test !occursin("`$(field)`", result_text)
        end
        @test occursin("ResolvedTarget, ENSG00000120948.20, ENST00000240185.8",
            result_text)
        @test occursin(
            "ThorAxeMSAResult, 1 seed, selected PID 10.0 (100 seqs, 1000 cols)",
            result_text)
        @test occursin(
            "1 expansion, 1 completed, 0 hits, selected MSA (250 seqs, 1000 cols)",
            result_text)
        @test occursin("1 result, ok", result_text)
        @test occursin("0 stages, empty", result_text)
        @test occursin("1 warning", result_text)
        @test occursin(":warn", result_text)
        @test !occursin("result.status", result_text)
        @test !occursin("result.target", result_text)
        @test !occursin("ResolvedTarget(", result_text)
        @test !occursin("ThorAxeMSAResult(", result_text)
        @test !occursin("ExpansionResult(", result_text)
        @test !occursin("ValidationResult(", result_text)
        @test !occursin("identity_hash", result_text)
        @test !occursin("\e[", result_text)

        missing_result_text = repr("text/plain", missing_expansion_result)
        @test occursin("1 expansion, 1 missing, 0 hits", missing_result_text)
        @test occursin("0 results, empty", missing_result_text)
        @test occursin(":error", missing_result_text)

        function _display_expansion(label; status = :ok)
            return Iduna.ExpansionResult(;
                run_dir = label,
                seed_stockholm = "seed.sto",
                hits_fasta = "hits.fasta",
                full_stockholm = "full.sto",
                match_stockholm = "match.sto",
                a3m_path = "$(label).a3m",
                db_dir = "db",
                hmm_dir = "hmm",
                logs_dir = "logs/expansion",
                status)
        end

        empty_target = Iduna.ResolvedTarget(;
            input_id = "empty",
            input_kind = :uniprot,
            ensembl_gene_id = "",
            transcript_id = "")
        empty_thoraxe = Iduna.ThorAxeMSAResult(;
            input_dir = "thoraxe_input",
            thoraxe_dirs = String[],
            msa_dir = "thoraxe_msa",
            baseline_fastas = String[],
            baseline_stockholms = String[],
            sequence_fastas = String[],
            species_files = String[],
            pid_summary = "candidate_summary.csv",
            seeds = Iduna.SeedSelection[],
            logs_dir = "logs/thoraxe")
        empty_sections_result = Iduna.IdunaResult(;
            input_id = "empty",
            workdir = "workdir",
            target = empty_target,
            thoraxe_msa = empty_thoraxe,
            expansions = Union{Missing, Iduna.ExpansionResult}[],
            validations = Iduna.ValidationResult[],
            stages = Any["no status"],
            status = :unknown)
        empty_sections_text = repr("text/plain", empty_sections_result)
        @test occursin("ResolvedTarget, unknown, unknown", empty_sections_text)
        @test occursin("ThorAxeMSAResult, 0 seeds, selected PIDs unknown",
            empty_sections_text)
        @test occursin("0 expansions, empty, 0 hits", empty_sections_text)
        @test occursin("1 stage, 1 unknown", empty_sections_text)

        warn_expansion = _display_expansion("warn_expansion"; status = :warn)
        warn_validation = Iduna.ValidationResult(;
            stats_path = "warn_stats.csv",
            status = :warn)
        warn_count_result = Iduna.IdunaResult(;
            input_id = "Q13148",
            workdir = "workdir",
            target,
            thoraxe_msa = thoraxe,
            expansions = [warn_expansion],
            validations = [warn_validation],
            stages = Any[Dict("status" => "warn")],
            status = :warn)
        warn_count_text = repr("text/plain", warn_count_result)
        @test occursin("1 expansion, 1 warn, 0 hits", warn_count_text)
        @test occursin("ThorAxeMSAResult, 1 seed, selected PID 10.0",
            warn_count_text)
        @test !occursin("(seqs, cols)", warn_count_text)
        @test occursin("1 result, warn", warn_count_text)
        @test occursin("1 stage, 1 warn", warn_count_text)

        second_seed = Iduna.SeedSelection(;
            pid = 80.0,
            epli = 0.8,
            stockholm_path = "seed80.sto",
            summary_path = "candidate_summary.csv")
        multi_thoraxe = Iduna.ThorAxeMSAResult(;
            input_dir = "thoraxe_input",
            thoraxe_dirs = ["thoraxe10", "thoraxe80"],
            msa_dir = "thoraxe_msa",
            baseline_fastas = ["msa10.fasta", "msa80.fasta"],
            baseline_stockholms = ["msa10.sto", "msa80.sto"],
            sequence_fastas = ["msa10_sequences.fasta", "msa80_sequences.fasta"],
            species_files = ["species10.txt", "species80.txt"],
            pid_summary = "candidate_summary.csv",
            seeds = [best_seed, second_seed],
            logs_dir = "logs/thoraxe")
        multi_result = Iduna.IdunaResult(;
            input_id = "Q13148",
            workdir = "workdir",
            target,
            thoraxe_msa = multi_thoraxe,
            expansions = [expansion, _display_expansion("expansion80")],
            validations = [
                validation,
                Iduna.ValidationResult(;
                    stats_path = "stats80.csv",
                    seed_nseq = 95,
                    seed_ncol = 900,
                    expanded_nseq = 275,
                    expanded_ncol = 900)
            ],
            status = :ok)
        multi_result_text = repr("text/plain", multi_result)
        @test occursin(
            "ThorAxeMSAResult, 2 seeds, selected PIDs (seqs, cols) 10.0 (100, 1000), 80.0 (95, 900)",
            multi_result_text)
        @test occursin(
            "2 expansions, 2 completed, 0 hits, selected MSAs (seqs, cols) 10.0 (250, 1000), 80.0 (275, 900)",
            multi_result_text)

        unknown_status_result = Iduna.IdunaResult(;
            input_id = "Q13148",
            workdir = "workdir",
            target,
            thoraxe_msa = thoraxe,
            expansions = Union{Missing, Iduna.ExpansionResult}[missing],
            validations = Iduna.ValidationResult[],
            stages = Any[(; status = :paused)],
            status = :mystery)
        unknown_result_text = repr("text/plain", unknown_status_result)
        @test occursin("1 unknown", unknown_result_text)
        @test occursin(":mystery", unknown_result_text)

        unknown_expansion = _display_expansion("unknown_expansion"; status = :mystery)
        unknown_expansion_result = Iduna.IdunaResult(;
            input_id = "Q13148",
            workdir = "workdir",
            target,
            thoraxe_msa = thoraxe,
            expansions = [unknown_expansion],
            validations = Iduna.ValidationResult[],
            status = :unknown)
        unknown_expansion_text = repr("text/plain", unknown_expansion_result)
        @test occursin("1 expansion, 1 unknown, 0 hits", unknown_expansion_text)

        function _colored_result_text(displayed_result)
            buffer = IOBuffer()
            show(IOContext(buffer, :color => true), MIME"text/plain"(), displayed_result)
            return String(take!(buffer))
        end

        colored_result_text = _colored_result_text(result)
        @test occursin("\e[1mFIELD", colored_result_text)
        @test occursin("FIELD      \e[22m", colored_result_text)
        @test occursin("\e[1mSUMMARY\e[22m", colored_result_text)
        @test occursin("\e[33m:warn\e[39m", colored_result_text)
        @test occursin("\e[33m1 warning\e[39m", colored_result_text)
        @test occursin("\e[32mok\e[39m", colored_result_text)
        @test occursin("\e[32m1 completed\e[39m", colored_result_text)

        done_result = Iduna.IdunaResult(;
            input_id = "Q13148",
            workdir = "workdir",
            target,
            thoraxe_msa = thoraxe,
            expansions = [expansion],
            validations = [validation],
            stages = Any[(; status = :done)],
            status = :ok)
        done_result_text = _colored_result_text(done_result)
        @test occursin("\e[32m:ok\e[39m", done_result_text)
        @test occursin("\e[32m1 done\e[39m", done_result_text)

        failed_expansion = _display_expansion("failed_expansion"; status = :failed)
        failed_validation = Iduna.ValidationResult(;
            stats_path = "failed_stats.csv",
            status = :failed)
        failed_result = Iduna.IdunaResult(;
            input_id = "Q13148",
            workdir = "workdir",
            target,
            thoraxe_msa = thoraxe,
            expansions = [failed_expansion],
            validations = [failed_validation],
            stages = Any[(; status = :failed)],
            status = :error)
        failed_result_text = _colored_result_text(failed_result)
        @test occursin("\e[31m:error\e[39m", failed_result_text)
        @test occursin("\e[31m1 failed\e[39m", failed_result_text)

        unknown_result_colored_text = _colored_result_text(unknown_status_result)
        @test occursin("\e[90m1 missing\e[39m", unknown_result_colored_text)
        @test occursin("\e[90mempty\e[39m", unknown_result_colored_text)
        @test occursin("\e[90m1 unknown\e[39m", unknown_result_colored_text)
        @test occursin("\e[90m:mystery\e[39m", unknown_result_colored_text)
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
                epli = 1.0,
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
            @test !haskey(written, "workdir")
            @test written["target"]["uniprot_sequence_path"] ==
                  joinpath("sequences", "uniprot", "Q13148.fasta")
            @test written["thoraxe_msa"]["pid_summary"] ==
                  joinpath("thoraxe_msa", "candidate_summary.csv")
            @test written["thoraxe_msa"]["seeds"][1]["stockholm_path"] ==
                  result.thoraxe_msa.seeds[1].stockholm_path
            @test written["thoraxe_msa"]["seeds"][1]["summary_path"] ==
                  result.thoraxe_msa.seeds[1].summary_path
            @test written["expansions"][1]["match_stockholm"] ==
                  result.expansions[1].match_stockholm
            @test written["expansions"][1]["s_exon_blocks_tsv"] ==
                  result.expansions[1].s_exon_blocks_tsv
            @test written["validations"][1]["stats_path"] ==
                  joinpath("validation", "pid_10.00", "stats.csv")
            result_state = JSON.parse(read(
                joinpath(workdir, ".iduna", "stages", "result", "stage_state.json"),
                String))
            @test result_state["outputs"]["result_json"] == "result.json"
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
                mmseqs_threads = 7,
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
            @test captured[][:mmseqs_threads] == 7
        end
    end

    @testset "generic threads keyword is not accepted" begin
        @test_throws MethodError Iduna.iduna(; id = "Q13148", threads = 1)
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
                sampling_strategy = :input,
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
            @test captured[][:sampling_strategy] === :input
        end
    end

    @testset "specieslist option forwarding" begin
        mktempdir() do tmp
            captured_default = Ref{Dict{Symbol, Any}}()
            Iduna.iduna(;
                id = "Q13148",
                workdir = joinpath(tmp, "default_specieslist_forwarding"),
                no_expansion = true,
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (
                    args...; kwargs...) -> begin
                    captured_default[] = Dict{Symbol, Any}(kwargs)
                    thoraxe
                end,
                _validate_results = (args...; kwargs...) -> validation)
            @test captured_default[][:specieslist] == "ases"

            captured_all = Ref{Dict{Symbol, Any}}()
            Iduna.iduna(;
                id = "Q13148",
                workdir = joinpath(tmp, "all_specieslist_forwarding"),
                no_expansion = true,
                specieslist = "all",
                _resolve_target = (args...; kwargs...) -> target,
                _build_thoraxe_msa = (
                    args...; kwargs...) -> begin
                    captured_all[] = Dict{Symbol, Any}(kwargs)
                    thoraxe
                end,
                _validate_results = (args...; kwargs...) -> validation)
            @test captured_all[][:specieslist] == "all"
        end
    end

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
end
