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
