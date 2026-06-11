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
