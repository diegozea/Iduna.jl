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
        summary_path = "best_seed.csv")
    thoraxe = Iduna.ThorAxeMSAResult(;
        input_dir = "thoraxe_input",
        thoraxe_dir = "thoraxe",
        msa_dir = "thoraxe_msa",
        baseline_fasta = "msa_0.fasta",
        baseline_stockholm = "msa_0.sto",
        sequence_fasta = "msa_0_sequences.fasta",
        species_file = "species.txt",
        pid_summary = "best_seed.csv",
        best_seed,
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
        expansion,
        validation,
        warnings = thoraxe.warnings,
        status = :warn)
    summary = Iduna.Utils.result_summary(result)
    @test summary.status == "warn"
    @test summary.thoraxe_msa.status == "warn"
    @test summary.thoraxe_msa.warnings == thoraxe.warnings
end
