mktempdir() do tmp
    fasta = joinpath(tmp, "candidate.fasta")
    write(fasta, ">ENSG\nA-C\n>ORTHO1\nABC\n")
    msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
    Iduna.Utils.set_s_exon_annotations!(msa, "000", ['0' => "1_0"])
    paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 60.0, 0)
    Iduna.ThorAxeMSA._write_candidate_sample_inputs(
        paths, msa, ["homo_sapiens", "mus_musculus"], [1, 2])
    Iduna.ThorAxeMSA.write_file(paths.fasta_path, msa, Iduna.ThorAxeMSA.FASTA)
    Iduna.ThorAxeMSA.write_file(
        paths.stockholm_path, msa, Iduna.ThorAxeMSA.Stockholm)
    written_species = read(paths.species_file, String)
    @test Iduna.ThorAxeMSA._write_species_file(
        paths.species_file, ["danio_rerio"]) == paths.species_file
    @test read(paths.species_file, String) == written_species
    @test Iduna.ThorAxeMSA._write_candidate_sample_inputs(
        paths, msa, ["wrong_length"], [1]) ==
          (paths.sequence_fasta, paths.species_file)
    @test_throws ErrorException Iduna.ThorAxeMSA._write_candidate_sample_inputs(
        Iduna.ThorAxeMSA._pid_sample_paths(tmp, 70.0, 0),
        msa, ["wrong_length"], [1])

    thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 60.0, 0)
    mkpath(thoraxe_dir)
    write(joinpath(thoraxe_dir, "path_table.csv"), "TranscriptIDCluster,Path\n")
    @test !Iduna.ThorAxeMSA._has_phylosofs_outputs(thoraxe_dir)
    phylosofs_dir = joinpath(thoraxe_dir, "phylosofs")
    mkpath(phylosofs_dir)
    write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\t0\n")
    write(joinpath(phylosofs_dir, "transcripts.pir"), ">P1;ENST\n0\nA*\n")
    @test Iduna.ThorAxeMSA._has_phylosofs_outputs(thoraxe_dir)
    @test Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
        Iduna.ResolvedTarget(;
            input_id = "ENST",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG",
            transcript_id = "ENST"),
        joinpath(tmp, "input"), tmp, 60.0, nothing, 0;
        keep_thoraxe_dir = true) ==
          (paths.fasta_path, paths.stockholm_path, thoraxe_dir)

    kept_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 90.0, 0)
    kept_thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 90.0, 0)
    mkpath(dirname(kept_paths.fasta_path))
    write_fake_thoraxe_dir(kept_thoraxe_dir)
    @test Iduna.ThorAxeMSA._run_kept_thoraxe_pid_msa!(
        Iduna.ResolvedTarget(;
            input_id = "ENST",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG",
            transcript_id = "ENST"),
        joinpath(tmp, "input"), tmp, kept_paths, kept_thoraxe_dir,
        joinpath(kept_thoraxe_dir, "path_table.csv"), 90.0, nothing, 0,
        command -> nothing, false) ==
          (kept_paths.fasta_path, kept_paths.stockholm_path, kept_thoraxe_dir)
    @test isfile(kept_paths.stockholm_path)

    fake_calls = Ref(0)
    fake_thoraxe = (input_dir, run_root; identity, specieslist, phylosofs,
        runner) -> begin
        fake_calls[] += 1
        @test identity == 50.0
        @test phylosofs === true
        write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
        nothing
    end
    temp_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 50.0, 0)
    pid_specieslist = "homo_sapiens,mus_musculus,danio_rerio," *
                      "xenopus_tropicalis,gallus_gallus,canis_lupus"
    pid_logs,
    pid_result = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
            Iduna.ResolvedTarget(;
                input_id = "ENST",
                input_kind = :ensembl_transcript,
                ensembl_gene_id = "ENSG",
                transcript_id = "ENST"),
            joinpath(tmp, "input"), tmp, 50.0, pid_specieslist, 0;
            thoraxe_fn = fake_thoraxe)
    end
    @test pid_result ==
          (temp_paths.fasta_path, temp_paths.stockholm_path,
        Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 50.0, 0))
    @test isempty([log
                   for log in pid_logs
                   if log.message in ("Running ThorAxe PID MSA.",
        "Running ThorAxe PID MSA command.")])
    @test fake_calls[] == 1
    @test isfile(temp_paths.stockholm_path)

    generated_pid = 55.0
    generated_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, generated_pid, 0)
    generated_thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(
        tmp, generated_pid, 0)
    write_fake_thoraxe_dir(generated_thoraxe_dir)
    generated_candidate = Iduna.ThorAxeMSA._generate_pid_candidate(
        Iduna.ResolvedTarget(;
            input_id = "ENST",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG",
            transcript_id = "ENST"),
        joinpath(tmp, "input"), tmp, generated_pid, nothing)
    @test generated_candidate.fasta_path == generated_paths.fasta_path
    @test generated_candidate.stockholm_path == generated_paths.stockholm_path
    @test generated_candidate.thoraxe_dir == generated_thoraxe_dir
    @test generated_candidate.s_exon_blocks_tsv == generated_paths.s_exon_blocks_tsv
    @test isfile(generated_candidate.stockholm_path)
    @test isfile(generated_candidate.s_exon_blocks_tsv)
    generated_sequence_fasta = read(generated_candidate.sequence_fasta, String)
    @test occursin(">ENSG\nAA", generated_sequence_fasta)
    @test occursin(">ORTHO1\nAX", generated_sequence_fasta)
    @test read(generated_candidate.species_file, String) ==
          "homo_sapiens\nmus_musculus\n"

    scoring_target = Iduna.ResolvedTarget(;
        input_id = "ENST",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG",
        transcript_id = "ENST")
    scoring_input = write_test_ensembl_bundle(joinpath(tmp, "scoring_input"))
    zero_sample_pid = 52.0
    write_fake_thoraxe_dir(
        Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, zero_sample_pid, 0))
    zero_sample_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        scoring_input, scoring_target, [zero_sample_pid];
        sample_count = 0,
        sample_fraction = 1.0,
        sample_seed = UInt64(11),
        requested_sample_seed = 11,
        sampling_strategy = :independent,
        effective_specieslist = nothing,
        orthology = "1:1",
        specieslist_filter = false,
        biomart_datasets_filter = false)
    zero_sample_row = Iduna.ThorAxeMSA._score_pid_candidate(
        scoring_target, scoring_input, tmp, zero_sample_pid, 1, nothing;
        sample_count = 0,
        sample_fraction = 1.0,
        sample_seed = UInt64(11),
        metadata = zero_sample_metadata)
    @test zero_sample_row.pid == zero_sample_pid
    @test zero_sample_row.eligible
    @test zero_sample_row.selection_mode == "all_candidates"
    @test zero_sample_row.n_samples == 0

    invalid_uniprot = joinpath(tmp, "invalid_uniprot.fasta")
    write(invalid_uniprot, ">P0\nA\n")
    invalid_target = Iduna.ResolvedTarget(;
        input_id = "P0",
        input_kind = :uniprot,
        uniprot_id = "P0",
        ensembl_gene_id = "ENSG",
        transcript_id = "ENST",
        uniprot_sequence_path = invalid_uniprot)
    invalid_pid = 52.5
    write_fake_thoraxe_dir(
        Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, invalid_pid, 0))
    invalid_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        scoring_input, invalid_target, [invalid_pid];
        sample_count = 0,
        sample_fraction = 1.0,
        sample_seed = UInt64(13),
        requested_sample_seed = 13,
        sampling_strategy = :independent,
        effective_specieslist = nothing,
        orthology = "1:1",
        specieslist_filter = false,
        biomart_datasets_filter = false)
    invalid_row = @test_logs (:info, r"Skipping ineligible ThorAxe PID candidate") match_mode=:any Iduna.ThorAxeMSA._score_pid_candidate(
        invalid_target, scoring_input, tmp, invalid_pid, 1, nothing;
        sample_count = 0,
        sample_fraction = 1.0,
        sample_seed = UInt64(13),
        metadata = invalid_metadata)
    @test invalid_row.pid == invalid_pid
    @test !invalid_row.eligible
    @test invalid_row.msa0_status == "invalid_msa0"
    @test occursin("indels", invalid_row.msa0_issue)

    invalid_two_step_pid = 52.75
    invalid_two_step_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        scoring_input, invalid_target, [invalid_two_step_pid];
        sample_count = 0,
        sample_fraction = 1.0,
        sample_seed = UInt64(14),
        requested_sample_seed = 14,
        sampling_strategy = :independent,
        effective_specieslist = nothing,
        orthology = "1:1",
        specieslist_filter = false,
        biomart_datasets_filter = false)
    invalid_two_step_calls = Ref(0)
    invalid_two_step_thoraxe = (input_dir, run_root; identity, specieslist,
        phylosofs, runner) -> begin
        invalid_two_step_calls[] += 1
        write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
        nothing
    end
    invalid_two_step_rows = @test_logs (
        :info, r"Skipping ineligible ThorAxe PID candidate") match_mode=:any Iduna.ThorAxeMSA._score_pid_candidates(
        invalid_target, scoring_input, tmp, [invalid_two_step_pid], nothing,
        invalid_two_step_metadata;
        pid_sample_count = 0,
        pid_sample_fraction = 1.0,
        sample_seed = UInt64(14),
        overwrite = true,
        thoraxe_fn = invalid_two_step_thoraxe)
    invalid_two_step_row = only(invalid_two_step_rows)
    @test invalid_two_step_calls[] == 1
    @test invalid_two_step_row.pid == invalid_two_step_pid
    @test !invalid_two_step_row.eligible
    @test invalid_two_step_row.msa0_status == "invalid_msa0"
    @test occursin("indels", invalid_two_step_row.msa0_issue)

    sampled_pid = 53.0
    sampled_calls = Ref(0)
    sampled_thoraxe = (input_dir, run_root; identity, specieslist, phylosofs,
        runner) -> begin
        sampled_calls[] += 1
        write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
        nothing
    end
    sampled_score = (reference_fasta, sample_fasta; logs_dir,
        label) -> begin
        @test isfile(reference_fasta)
        @test isfile(sample_fasta)
        @test label == "species_subset_001"
        mkpath(logs_dir)
        write(joinpath(logs_dir, "$(label)_hhalign.out"), "fake\n")
        (; raw_score = 88.0, matched_positions = 88, comparable_positions = 100)
    end
    sampled_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        scoring_input, scoring_target, [sampled_pid];
        sample_count = 1,
        sample_fraction = 1.0,
        sample_seed = UInt64(12),
        requested_sample_seed = 12,
        sampling_strategy = :independent,
        effective_specieslist = nothing,
        orthology = "1:1",
        specieslist_filter = false,
        biomart_datasets_filter = false)
    sampled_progress_stderr = joinpath(tmp, "sampled_progress_stderr.log")
    sampled_logs,
    sampled_row = Test.collect_test_logs() do
        open(sampled_progress_stderr, "w") do stderr_io
            redirect_stderr(stderr_io) do
                Iduna.ThorAxeMSA._score_pid_candidate(
                    scoring_target, scoring_input, tmp, sampled_pid, 1, nothing;
                    sample_count = 1,
                    sample_fraction = 1.0,
                    sample_seed = UInt64(12),
                    metadata = sampled_metadata,
                    thoraxe_fn = sampled_thoraxe,
                    score_fn = sampled_score)
            end
        end
    end
    @test !any(
        log -> log.message in ("Scoring ThorAxe PID sample.",
            "Scoring ThorAxe PID candidate.",
            "Preparing ThorAxe PID candidate."),
        sampled_logs)
    @test isempty(read(sampled_progress_stderr, String))
    @test sampled_calls[] == 2
    @test sampled_row.epli == 88.0
    @test sampled_row.n_samples == 1
    @test isfile(Iduna.ThorAxeMSA._pid_scores_path(tmp, sampled_pid))

    progress_pid = 53.25
    progress_ref_fasta = joinpath(tmp, "progress_reference.fasta")
    write(progress_ref_fasta, ">homo_sapiens\nAA\n>mus_musculus\nAA\n")
    progress_candidate = (;
        msa = Iduna.ThorAxeMSA.read_file(
            progress_ref_fasta, Iduna.ThorAxeMSA.FASTA),
        fasta_path = progress_ref_fasta)
    progress_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, progress_pid, 1)
    mkpath(dirname(progress_paths.species_file))
    write(progress_paths.species_file, "homo_sapiens\n")
    progress_thoraxe = (input_dir, run_root; identity, specieslist, phylosofs,
        runner) -> begin
        write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
        nothing
    end
    progress_output = IOBuffer()
    progress_logs,
    progress_rows = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._score_pid_samples(
            scoring_target, scoring_input, tmp, progress_candidate,
            progress_pid, 1;
            sample_count = 2,
            thoraxe_fn = progress_thoraxe,
            score_fn = sampled_score,
            progress_output,
            progress_enabled = true)
    end
    progress_text = String(take!(progress_output))
    @test occursin("  0%", progress_text)
    @test occursin("100%", progress_text)
    @test !occursin("Scoring ThorAxe PID samples in parallel", progress_text)
    @test !any(log -> log.message == "Scoring ThorAxe PID samples in parallel.",
        progress_logs)
    @test length(progress_rows.rows) == 1

    parallel_pid = 53.5
    parallel_lock = ReentrantLock()
    parallel_write_lock = ReentrantLock()
    parallel_call_kinds = Symbol[]
    parallel_sample_threads = Int[]
    parallel_score_labels = String[]
    function record_parallel_call!(kind::Symbol)
        lock(parallel_lock)
        try
            push!(parallel_call_kinds, kind)
            if kind === :sample
                push!(parallel_sample_threads, Base.Threads.threadid())
            end
        finally
            unlock(parallel_lock)
        end
    end
    parallel_thoraxe = (input_dir, run_root; identity, specieslist, phylosofs,
        runner) -> begin
        record_parallel_call!(specieslist === nothing ? :full : :sample)
        specieslist === nothing || sleep(0.02)
        lock(parallel_write_lock)
        try
            write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
        finally
            unlock(parallel_write_lock)
        end
        nothing
    end
    parallel_score = (reference_fasta, sample_fasta; logs_dir,
        label) -> begin
        lock(parallel_lock)
        try
            push!(parallel_score_labels, String(label))
        finally
            unlock(parallel_lock)
        end
        mkpath(logs_dir)
        write(joinpath(logs_dir, "$(label)_hhalign.out"), "fake\n")
        value = parse(Float64, last(split(String(label), '_')))
        (; raw_score = value, matched_positions = value, comparable_positions = 100)
    end
    parallel_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        scoring_input, scoring_target, [parallel_pid];
        sample_count = 4,
        sample_fraction = 1.0,
        sample_seed = UInt64(16),
        requested_sample_seed = 16,
        sampling_strategy = :independent,
        effective_specieslist = nothing,
        orthology = "1:1",
        specieslist_filter = false,
        biomart_datasets_filter = false)
    parallel_row = Iduna.ThorAxeMSA._score_pid_candidate(
        scoring_target, scoring_input, tmp, parallel_pid, 1, nothing;
        sample_count = 4,
        sample_fraction = 1.0,
        sample_seed = UInt64(16),
        metadata = parallel_metadata,
        thoraxe_fn = parallel_thoraxe,
        score_fn = parallel_score)
    @test count(==(:full), parallel_call_kinds) == 1
    @test count(==(:sample), parallel_call_kinds) == 4
    @test sort(parallel_score_labels) == ["species_subset_001",
        "species_subset_002", "species_subset_003", "species_subset_004"]
    if Base.Threads.threadpoolsize() > 1
        @test length(unique(parallel_sample_threads)) > 1
    end
    @test parallel_row.epli == 2.5
    @test parallel_row.n_samples == 4
    parallel_scores = Iduna.ThorAxeMSA.DataFrame(
        Iduna.ThorAxeMSA.CSV.File(
        Iduna.ThorAxeMSA._pid_scores_path(tmp, parallel_pid)))
    @test parallel_scores.sample == 1:4

    common_pids = [54.0, 56.0]
    common_calls = Ref(0)
    common_thoraxe = (input_dir, run_root; identity, specieslist, phylosofs,
        runner) -> begin
        common_calls[] += 1
        write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
        nothing
    end
    common_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
        scoring_input, scoring_target, common_pids;
        sample_count = 1,
        sample_fraction = 1.0,
        sample_seed = UInt64(15),
        requested_sample_seed = 15,
        sampling_strategy = :common,
        effective_specieslist = nothing,
        orthology = "1:1",
        specieslist_filter = false,
        biomart_datasets_filter = false)
    common_logs,
    common_rows = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._score_pid_candidates(
            scoring_target, scoring_input, tmp, common_pids, nothing,
            common_metadata;
            pid_sample_count = 1,
            pid_sample_fraction = 1.0,
            sample_seed = UInt64(15),
            overwrite = true,
            thoraxe_fn = common_thoraxe,
            score_fn = sampled_score)
    end
    @test isempty([log
                   for log in common_logs
                   if log.message == "Preparing ThorAxe PID candidate."])
    common_sampling_log = only([log
                                for log in common_logs
                                if log.message ==
                                   "Preparing common ThorAxe PID species samples."])
    common_sampling_kwargs = Dict(common_sampling_log.kwargs)
    @test common_sampling_kwargs[:n_common_species] == 2
    @test !haskey(common_sampling_kwargs, :sample_count)
    @test common_calls[] == 4
    @test [row.pid for row in common_rows] == common_pids
    @test all(row -> row.n_samples == 1, common_rows)
    common_species_file = Iduna.ThorAxeMSA._shared_sample_species_file(tmp, 1)
    @test isfile(common_species_file)
    for pid in common_pids
        species_file = Iduna.ThorAxeMSA._pid_sample_paths(tmp, pid, 1).species_file
        @test islink(species_file)
        @test read(species_file, String) == read(common_species_file, String)
    end

    scored_rows = Iduna.ThorAxeMSA._score_pid_candidates(
        scoring_target, scoring_input, tmp, [62.0, 63.0], nothing,
        zero_sample_metadata;
        pid_sample_count = 0,
        pid_sample_fraction = 1.0,
        sample_seed = UInt64(13),
        overwrite = false,
        thoraxe_fn = sampled_thoraxe)
    @test [row.pid for row in scored_rows] == [62.0, 63.0]
    @test all(row -> row.selection_mode == "all_candidates", scored_rows)

    stale_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 80.0, 0)
    mkpath(dirname(stale_paths.fasta_path))
    Iduna.ThorAxeMSA.write_file(stale_paths.fasta_path, msa, Iduna.ThorAxeMSA.FASTA)
    Iduna.ThorAxeMSA.write_file(
        stale_paths.stockholm_path, msa, Iduna.ThorAxeMSA.Stockholm)
    stale_thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 80.0, 0)
    mkpath(stale_thoraxe_dir)
    write(joinpath(stale_thoraxe_dir, "path_table.csv"), "TranscriptIDCluster,Path\n")
    @test_throws Exception Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
        Iduna.ResolvedTarget(;
            input_id = "ENST",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG",
            transcript_id = "ENST"),
        joinpath(tmp, "missing_input"), tmp, 80.0, nothing, 0;
        keep_thoraxe_dir = true)
    @test_throws Exception Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
        Iduna.ResolvedTarget(;
            input_id = "ENST",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG",
            transcript_id = "ENST"),
        joinpath(tmp, "missing_input"), tmp, 60.0, nothing, 0;
        overwrite = true,
        keep_thoraxe_dir = true)

    positions,
    codes = Iduna.ThorAxeMSA._get_codes(
        "Probab=99\nquery 1 A-C 2\n        | |\nConfidence\n")
    @test positions == [1, 0, 2]
    @test codes == ['|', ' ', '|']
    @test Iduna.ThorAxeMSA._identity_from_codes(positions, codes) == 100.0
    @test_throws ErrorException Iduna.ThorAxeMSA._reference_index(
        msa, "MISSING_GENE", "MISSING_TRANSCRIPT")

    reference_fasta = joinpath(tmp, "identity_reference.fasta")
    sample_fasta = joinpath(tmp, "identity_sample.fasta")
    write(reference_fasta, ">ref\nAA\n")
    write(sample_fasta, ">sample\nAA\n")
    identity_logs = joinpath(tmp, "identity_logs")
    @test Iduna.ThorAxeMSA.compute_identity_against_reference(
        reference_fasta, sample_fasta; logs_dir = identity_logs,
        label = "identical") == 100.0
    @test isfile(joinpath(identity_logs, "identical_hhalign.out"))
end
