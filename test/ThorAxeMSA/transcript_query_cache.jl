mktempdir() do tmp
    target = Iduna.ResolvedTarget(;
        input_id = "ENST",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG",
        transcript_id = "ENST",
        species = "homo_sapiens")
    fasta = joinpath(tmp, "candidate.fasta")
    write(fasta,
        ">ENSG\nAAAA\n" *
        ">ORTHO1\nBBBB\n")
    msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
    records = [
        (;
            pid = 10.0,
            candidate = (;
                msa,
                species = ["homo_sapiens", "mus_musculus"]),
            validation = (; eligible = true)),
        (;
            pid = 20.0,
            candidate = (;
                msa,
                species = ["homo_sapiens", "mus_musculus"]),
            validation = (; eligible = true))
    ]
    input_logs,
    _ = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._prepare_candidate_species_samples!(
            target, joinpath(tmp, "input"), tmp, records,
            "homo_sapiens,mus_musculus,danio_rerio", :input;
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = UInt64(9),
            overwrite = true)
    end
    @test isempty([log
                   for log in input_logs
                   if log.message == "Preparing common ThorAxe PID species samples."])
    input_info = only([log
                       for log in input_logs
                       if log.message == "Preparing input ThorAxe PID species samples."])
    input_info_kwargs = Dict(input_info.kwargs)
    @test input_info_kwargs[:n_species] == 3
    @test !haskey(input_info_kwargs, :sample_count)
    canonical = Iduna.ThorAxeMSA._shared_sample_species_file(tmp, 1)
    @test Set(split(chomp(read(canonical, String)), '\n')) ==
          Set(["homo_sapiens", "mus_musculus", "danio_rerio"])
    for pid in (10.0, 20.0)
        link = Iduna.ThorAxeMSA._pid_sample_paths(tmp, pid, 1).species_file
        @test islink(link)
        @test read(link, String) == read(canonical, String)
    end
    @test_throws ErrorException Iduna.ThorAxeMSA._prepare_candidate_species_samples!(
        target, joinpath(tmp, "input"), tmp, records, nothing, :input;
        sample_count = 1,
        sample_fraction = 1.0,
        sample_seed = UInt64(9),
        overwrite = true)
end

mktempdir() do tmp
    input_dir = joinpath(tmp, "thoraxe_input")
    ensembl = joinpath(input_dir, "Ensembl")
    mkpath(ensembl)
    for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
        write(joinpath(ensembl, file), "$(file)\n")
    end
    metadata_target = Iduna.ResolvedTarget(;
        input_id = "P20963",
        input_kind = :uniprot,
        ensembl_gene_id = "ENSG",
        transcript_id = "ENST")
    metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        input_dir, metadata_target, [10.0, 80.0];
        sample_count = 2,
        sample_fraction = 0.5,
        sample_seed = UInt64(7),
        requested_sample_seed = 7,
        sampling_strategy = :independent,
        effective_specieslist = "homo_sapiens",
        orthology = "1:1",
        specieslist_filter = true,
        biomart_datasets_filter = false)
    summary = joinpath(tmp, "candidate_summary.csv")
    rows = [
        (;
            gene_id = "ENSG",
            transcript_id = "ENST",
            pid = 10.0,
            pid_order = 1,
            eligible = true,
            selected = true,
            msa0_status = "ok",
            msa0_issue = missing,
            epli = 70.0,
            n_samples = 2,
            n_sequences_msa0 = 4,
            pid_sample_count = metadata.pid_sample_count,
            pid_sample_fraction = metadata.pid_sample_fraction,
            pid_sample_seed = metadata.pid_sample_seed,
            sampling_strategy = String(metadata.sampling_strategy),
            pid_thresholds_key = metadata.pid_thresholds_key,
            effective_specieslist = metadata.effective_specieslist,
            orthology = metadata.orthology,
            specieslist_filter = metadata.specieslist_filter,
            biomart_datasets_filter = metadata.biomart_datasets_filter,
            transcript_query_fingerprint = metadata.transcript_query_fingerprint,
            selection_mode = metadata.selection_mode,
            fasta_path = "pid10.fasta",
            stockholm_path = "pid10.sto",
            sequence_fasta = "pid10_sequences.fasta",
            species_file = "pid10_species.txt",
            scores_path = "pid10_scores.csv"
        ),
        (;
            gene_id = "ENSG",
            transcript_id = "ENST",
            pid = 80.0,
            pid_order = 2,
            eligible = false,
            selected = false,
            msa0_status = "invalid_msa0",
            msa0_issue = "indels",
            epli = missing,
            n_samples = 0,
            n_sequences_msa0 = 3,
            pid_sample_count = metadata.pid_sample_count,
            pid_sample_fraction = metadata.pid_sample_fraction,
            pid_sample_seed = metadata.pid_sample_seed,
            sampling_strategy = String(metadata.sampling_strategy),
            pid_thresholds_key = metadata.pid_thresholds_key,
            effective_specieslist = metadata.effective_specieslist,
            orthology = metadata.orthology,
            specieslist_filter = metadata.specieslist_filter,
            biomart_datasets_filter = metadata.biomart_datasets_filter,
            transcript_query_fingerprint = metadata.transcript_query_fingerprint,
            selection_mode = metadata.selection_mode,
            fasta_path = "pid80.fasta",
            stockholm_path = "pid80.sto",
            sequence_fasta = "pid80_sequences.fasta",
            species_file = "pid80_species.txt",
            scores_path = "pid80_scores.csv"
        )
    ]
    Iduna.ThorAxeMSA.CSV.write(summary, Iduna.ThorAxeMSA.DataFrame(rows))
    df = Iduna.ThorAxeMSA._candidate_summary_dataframe(summary)
    @test Iduna.ThorAxeMSA._candidate_summary_matches(df, metadata)
    @test Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, metadata)
    @test Iduna.ThorAxeMSA._summary_sampling_strategy(df, :common) === :independent
    @test Iduna.ThorAxeMSA._summary_sampling_strategy(
        Iduna.ThorAxeMSA.DataFrame(pid = [10.0]), :common) === :common
    @test Iduna.ThorAxeMSA._summary_sampling_strategy(
        Iduna.ThorAxeMSA.DataFrame(sampling_strategy = String[]),
        :common) === :common
    @test Iduna.ThorAxeMSA._summary_sampling_strategy(
        Iduna.ThorAxeMSA.DataFrame(sampling_strategy = [missing]),
        :common) === :common
    changed_pids = merge(metadata, (; pid_thresholds_key = "80.0"))
    @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_pids)
    @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_pids)
    changed_fraction = merge(metadata, (; pid_sample_fraction = 0.75))
    @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_fraction)
    @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_fraction)
    changed_strategy = merge(metadata, (; sampling_strategy = :common))
    @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_strategy)
    @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_strategy)
    changed_gene = merge(metadata, (; gene_id = "ENSG_DIFFERENT"))
    @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_gene)
    @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_gene)
    changed_transcript = merge(metadata, (; transcript_id = "ENST_DIFFERENT"))
    @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_transcript)
    @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_transcript)
    @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(
        joinpath(tmp, "missing_summary.csv"), metadata)
    random_requested_seed = merge(metadata,
        (; pid_sample_seed = UInt64(99), requested_pid_sample_seed = nothing))
    @test Iduna.ThorAxeMSA._candidate_summary_matches(df, random_requested_seed)
    @test Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, random_requested_seed)
    @test Iduna.ThorAxeMSA._thoraxe_msa_identity(metadata, tmp) !=
          Iduna.ThorAxeMSA._thoraxe_msa_identity(
        merge(metadata,
            (; pid_sample_seed = UInt64(99), requested_pid_sample_seed = 99)),
        tmp)
    @test Iduna.ThorAxeMSA._thoraxe_msa_identity(random_requested_seed, tmp) ==
          Iduna.ThorAxeMSA._thoraxe_msa_identity(
        merge(random_requested_seed, (; pid_sample_seed = UInt64(7))), tmp)
    random_identity = Iduna.ThorAxeMSA._thoraxe_msa_identity(
        random_requested_seed, tmp)
    Iduna.ThorAxeMSA._write_thoraxe_msa_state(
        tmp, summary, Iduna.SeedSelection[], :done, random_identity; action = :run)
    random_stage_cache = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(
        tmp, summary, random_requested_seed,
        Iduna.ThorAxeMSA._thoraxe_msa_identity(
            merge(random_requested_seed, (; pid_sample_seed = UInt64(7))), tmp);
        overwrite = false)
    @test random_stage_cache.cache.reusable === true

    no_species_metadata = merge(metadata, (; effective_specieslist = nothing))
    fasta = joinpath(tmp, "candidate_summary_candidate.fasta")
    write(fasta,
        ">ENSG\nAAAA\n" *
        ">ORTHO1\nBBBB\n")
    candidate_msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
    candidate = (;
        msa = candidate_msa,
        fasta_path = "pid10.fasta",
        stockholm_path = "pid10.sto",
        sequence_fasta = "pid10_sequences.fasta",
        species_file = "pid10_species.txt",
        workdir = tmp)
    validation = (; eligible = true, status = "ok", issue = missing)
    no_species_rows = [
        Iduna.ThorAxeMSA._candidate_summary_row(
            metadata_target, candidate, 10.0, 1, validation;
            sample_count = no_species_metadata.pid_sample_count,
            sample_fraction = no_species_metadata.pid_sample_fraction,
            sample_seed = no_species_metadata.pid_sample_seed,
            metadata = no_species_metadata),
        Iduna.ThorAxeMSA._candidate_summary_row(
            metadata_target, candidate, 80.0, 2, validation;
            sample_count = no_species_metadata.pid_sample_count,
            sample_fraction = no_species_metadata.pid_sample_fraction,
            sample_seed = no_species_metadata.pid_sample_seed,
            metadata = no_species_metadata)
    ]
    no_species_summary = joinpath(tmp, "candidate_summary_no_species.csv")
    Iduna.ThorAxeMSA._summarize_candidate_scores(
        no_species_rows, no_species_summary)
    no_species_df = Iduna.ThorAxeMSA._candidate_summary_dataframe(no_species_summary)
    @test all(ismissing, no_species_df.effective_specieslist)
    @test Iduna.ThorAxeMSA._candidate_summary_matches(
        no_species_df, no_species_metadata)

    decimal_pid_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        input_dir, metadata_target, [12.34];
        sample_count = 2,
        sample_fraction = 0.5,
        sample_seed = UInt64(7),
        requested_sample_seed = 7,
        sampling_strategy = :independent,
        effective_specieslist = "homo_sapiens",
        orthology = "1:1",
        specieslist_filter = true,
        biomart_datasets_filter = false)
    @test decimal_pid_metadata.pid_thresholds_key == "12.34"
    decimal_pid_rows = [
        Iduna.ThorAxeMSA._candidate_summary_row(
        metadata_target, candidate, 12.34, 1, validation;
        epli = 70.0,
        n_samples = 2,
        sample_count = decimal_pid_metadata.pid_sample_count,
        sample_fraction = decimal_pid_metadata.pid_sample_fraction,
        sample_seed = decimal_pid_metadata.pid_sample_seed,
        metadata = decimal_pid_metadata)
    ]
    decimal_pid_summary = joinpath(tmp, "candidate_summary_decimal_pid.csv")
    Iduna.ThorAxeMSA._summarize_candidate_scores(
        decimal_pid_rows, decimal_pid_summary)
    decimal_pid_df = Iduna.ThorAxeMSA._candidate_summary_dataframe(
        decimal_pid_summary)
    @test Iduna.ThorAxeMSA._candidate_summary_matches(
        decimal_pid_df, decimal_pid_metadata)

    selected = Iduna.ThorAxeMSA._selected_candidate_seeds(summary)
    @test length(selected) == 1
    @test only(selected).pid == 10.0
    eligible = Iduna.ThorAxeMSA._eligible_candidate_seeds(summary)
    @test length(eligible) == 1
    @test only(eligible).pid == 10.0

    warning_summary = joinpath(tmp, "candidate_summary_warnings.csv")
    warning_rows = [
        Iduna.ThorAxeMSA._candidate_summary_row(
            metadata_target, candidate, 20.0, 1,
            (; eligible = false, status = "invalid_msa0", issue = "bad indel");
            sample_count = metadata.pid_sample_count,
            sample_fraction = metadata.pid_sample_fraction,
            sample_seed = metadata.pid_sample_seed,
            metadata),
        Iduna.ThorAxeMSA._candidate_summary_row(
            metadata_target, candidate, 30.0, 2,
            (; eligible = true, status = "warning", issue = "one substitution");
            sample_count = metadata.pid_sample_count,
            sample_fraction = metadata.pid_sample_fraction,
            sample_seed = metadata.pid_sample_seed,
            metadata)
    ]
    Iduna.ThorAxeMSA._summarize_candidate_scores(warning_rows, warning_summary)
    warning_text = Iduna.ThorAxeMSA._candidate_summary_warnings(warning_summary)
    @test any(w -> occursin("excluded", w) && occursin("bad indel", w),
        warning_text)
    @test any(
        w -> occursin("candidate retained", w) &&
             occursin("one substitution", w), warning_text)
end
