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
