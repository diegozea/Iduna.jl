mktempdir() do tmp
    workdir = joinpath(tmp, "iduna")
    summary = joinpath(workdir, "thoraxe_msa", "candidate_summary.csv")
    seed_rel = joinpath("thoraxe_msa", "candidates", "pid_10.00",
        "candidate_msa_full.sto")
    seed_fasta_rel = joinpath("thoraxe_msa", "candidates", "pid_10.00",
        "candidate_msa_full.fasta")
    seed_path = joinpath(workdir, seed_rel)
    seed_fasta = joinpath(workdir, seed_fasta_rel)
    mkpath(dirname(summary))
    mkpath(dirname(seed_path))
    write(seed_path, "# STOCKHOLM 1.0\nseed AC\n//\n")
    write(seed_fasta, ">seed\nAC\n")
    write(summary,
        "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
        "10.0,1,true,70.0,1,$(seed_rel),$(seed_fasta_rel)\n")

    seed = Iduna.ThorAxeMSA.select_best_seed(summary)
    @test seed.workdir == abspath(workdir)
    @test seed.stockholm_path == seed_rel
    outside = joinpath(tmp, "outside")
    mkpath(outside)
    cd(outside) do
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_seed_msa(seed)) == 1
    end
end

mktempdir() do tmp
    fasta = joinpath(tmp, "candidate.fasta")
    write(fasta,
        ">ORTHO1\nAAAA\n" *
        ">ENSG00000000001\nBBBB\n" *
        ">ORTHO2\nCCCC\n" *
        ">ORTHO3\nDDDD\n")
    msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
    species = ["mus_musculus", "homo_sapiens", "rattus_norvegicus", "danio_rerio"]

    rng = Iduna.ThorAxeMSA._sample_rng(UInt64(7), 1)
    indices = Iduna.ThorAxeMSA._sample_indices(4, 2, 0.5, rng)
    @test first(indices) == 2
    @test length(indices) == 3
    @test length(unique(indices)) == length(indices)
    @test !(2 in indices[2:end])

    full_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 80.0, 0)
    Iduna.ThorAxeMSA._write_candidate_sample_inputs(
        full_paths, msa, species, collect(1:4); overwrite = true)
    Iduna.ThorAxeMSA._ensure_pid_candidate_samples(tmp, 80.0, msa, species;
        sample_count = 1,
        sample_fraction = 0.5,
        sample_seed = UInt64(7),
        overwrite = true,
        gene_id = "ENSG00000000001.1",
        transcript_id = "ENST1")

    sample_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 80.0, 1)
    @test full_paths.species_file ==
          joinpath(tmp, "thoraxe_msa", "candidates", "pid_80.00", "species",
        "candidate_species_full.txt")
    @test sample_paths.species_file ==
          joinpath(tmp, "thoraxe_msa", "candidates", "pid_80.00", "species",
        "candidate_species_subset_001.txt")
    sample0 = split(chomp(read(full_paths.species_file, String)), '\n')
    sample1 = split(chomp(read(sample_paths.species_file, String)), '\n')
    @test sample0 == species
    @test first(sample1) == "homo_sapiens"
    @test length(sample1) == 3
    @test length(unique(sample1)) == length(sample1)
    @test all(item -> item in species[[1, 3, 4]], sample1[2:end])
end

mktempdir() do tmp
    target_path = joinpath(tmp, "shared", "candidate_species_subset_001.txt")
    mkpath(dirname(target_path))
    write(target_path, "homo_sapiens\n")
    link_path = joinpath(tmp, "pid", "species", "candidate_species_subset_001.txt")

    Iduna.ThorAxeMSA._symlink_species_file(link_path, target_path)
    @test islink(link_path)
    @test read(link_path, String) == "homo_sapiens\n"
    @test Iduna.ThorAxeMSA._symlink_species_file(link_path, target_path) ==
          link_path
    @test islink(link_path)

    regular_path = joinpath(tmp, "pid", "species", "existing_species.txt")
    mkpath(dirname(regular_path))
    write(regular_path, "existing\n")
    @test Iduna.ThorAxeMSA._symlink_species_file(regular_path, target_path) ==
          regular_path
    @test !islink(regular_path)
    @test read(regular_path, String) == "existing\n"
    Iduna.ThorAxeMSA._symlink_species_file(
        regular_path, target_path; overwrite = true)
    @test islink(regular_path)
    @test read(regular_path, String) == "homo_sapiens\n"

    existing_dir = joinpath(tmp, "pid", "species", "existing_dir")
    mkpath(existing_dir)
    @test Iduna.ThorAxeMSA._symlink_species_file(existing_dir, target_path) ==
          existing_dir
    @test isdir(existing_dir)

    rng = Iduna.ThorAxeMSA._sample_rng(UInt64(7), 1)
    @test Iduna.ThorAxeMSA._species_sample(
        ["homo_sapiens"], "homo_sapiens", 0.5, rng) == ["homo_sapiens"]
    @test_throws ErrorException Iduna.ThorAxeMSA._species_sample(
        ["mus_musculus"], "homo_sapiens", 1.0, rng)
end

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
    records_without_reference_overlap = [
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
                species = ["mus_musculus", "danio_rerio"]),
            validation = (; eligible = true))
    ]

    empty_common = Iduna.ThorAxeMSA._common_sampling_universe(
        NamedTuple[], "ENSG", "ENST")
    @test isempty(empty_common.species)
    @test empty_common.reference_species === nothing
    @test_throws ErrorException Iduna.ThorAxeMSA._common_sampling_universe(
        records_without_reference_overlap, "ENSG", "ENST")

    independent_workdir = joinpath(tmp, "independent")
    independent_logs,
    _ = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._prepare_candidate_species_samples!(
            target, joinpath(tmp, "input"), independent_workdir,
            [first(records_without_reference_overlap)], nothing, :independent;
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = UInt64(7),
            overwrite = true)
    end
    @test isempty([log
                   for log in independent_logs
                   if log.message == "Preparing common ThorAxe PID species samples."])
    independent_species_file = Iduna.ThorAxeMSA._pid_sample_paths(
        independent_workdir, 10.0, 1).species_file
    @test isfile(independent_species_file)
    @test !islink(independent_species_file)
    @test split(chomp(read(independent_species_file, String)), '\n') ==
          ["homo_sapiens", "mus_musculus"]

    input_universe = Iduna.ThorAxeMSA._input_sampling_universe(
        target, NamedTuple[], "homo_sapiens,mus_musculus")
    @test input_universe.species == ["homo_sapiens", "mus_musculus"]
    @test input_universe.reference_species == "homo_sapiens"
    @test_throws ErrorException Iduna.ThorAxeMSA._input_sampling_universe(
        target, NamedTuple[], "mus_musculus")
    @test_throws ErrorException Iduna.ThorAxeMSA._shared_sampling_universe(
        :independent, target, records_without_reference_overlap, nothing)

    for sampling_strategy in (:independent, :common, :input)
        @test Iduna.ThorAxeMSA._validate_sampling_strategy(sampling_strategy) ===
              sampling_strategy
    end
    @test_throws ErrorException Iduna.ThorAxeMSA._validate_sampling_strategy(:pid)
end

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
        ">ORTHO1\nBBBB\n" *
        ">ORTHO2\nCCCC\n" *
        ">ORTHO3\nDDDD\n")
    msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
    eligible = (; eligible = true)
    invalid = (; eligible = false)
    common_records = [
        (;
            pid = 10.0,
            candidate = (;
                msa,
                species = ["homo_sapiens", "mus_musculus", "danio_rerio",
                    "xenopus_tropicalis"]),
            validation = eligible),
        (;
            pid = 20.0,
            candidate = (;
                msa,
                species = ["homo_sapiens", "mus_musculus", "rattus_norvegicus",
                    "xenopus_tropicalis"]),
            validation = eligible),
        (;
            pid = 30.0,
            candidate = (;
                msa,
                species = ["homo_sapiens", "mus_musculus"]),
            validation = invalid)
    ]
    common_logs,
    _ = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._prepare_candidate_species_samples!(
            target, joinpath(tmp, "input"), tmp, common_records, nothing, :common;
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = UInt64(7),
            overwrite = true)
    end
    common_info = only([log
                        for log in common_logs
                        if log.message ==
                           "Preparing common ThorAxe PID species samples."])
    common_info_kwargs = Dict(common_info.kwargs)
    @test common_info_kwargs[:n_common_species] == 3
    @test !haskey(common_info_kwargs, :sample_count)
    common_warning = only([log
                           for log in common_logs
                           if log.message ==
                              "Common ThorAxe PID species list is small; consider sampling_strategy=:input or :independent."])
    common_warning_kwargs = Dict(common_warning.kwargs)
    @test common_warning_kwargs[:n_common_species] == 3
    @test common_warning_kwargs[:warning_threshold] == 6
    @test !haskey(common_warning_kwargs, :sample_count)
    canonical = Iduna.ThorAxeMSA._shared_sample_species_file(tmp, 1)
    common_species = split(chomp(read(canonical, String)), '\n')
    @test first(common_species) == "homo_sapiens"
    @test Set(common_species) ==
          Set(["homo_sapiens", "mus_musculus", "xenopus_tropicalis"])
    for pid in (10.0, 20.0)
        link = Iduna.ThorAxeMSA._pid_sample_paths(tmp, pid, 1).species_file
        @test islink(link)
        @test read(link, String) == read(canonical, String)
    end
    @test !ispath(Iduna.ThorAxeMSA._pid_sample_paths(tmp, 30.0, 1).species_file)

    provided_logs,
    _ = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._prepare_candidate_species_samples!(
            target, joinpath(tmp, "input"), joinpath(tmp, "provided_common"),
            common_records, "homo_sapiens,mus_musculus", :common;
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = UInt64(8),
            overwrite = true)
    end
    provided_warning = only([log
                             for log in provided_logs
                             if occursin("provided specieslist",
        String(log.message))])
    @test Dict(provided_warning.kwargs)[:n_common_species] == 3

    wide_species = ["homo_sapiens", "mus_musculus", "danio_rerio",
        "xenopus_tropicalis", "gallus_gallus", "canis_lupus", "bos_taurus"]
    wide_records = [
        (;
            pid = 40.0,
            candidate = (; msa, species = wide_species),
            validation = eligible),
        (;
            pid = 50.0,
            candidate = (; msa, species = wide_species),
            validation = eligible)
    ]
    wide_logs,
    _ = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._prepare_candidate_species_samples!(
            target, joinpath(tmp, "input"), joinpath(tmp, "wide_common"),
            wide_records, nothing, :common;
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = UInt64(9),
            overwrite = true)
    end
    wide_info = only([log
                      for log in wide_logs
                      if log.message ==
                         "Preparing common ThorAxe PID species samples."])
    wide_info_kwargs = Dict(wide_info.kwargs)
    @test wide_info_kwargs[:n_common_species] == 7
    @test !haskey(wide_info_kwargs, :sample_count)
    @test isempty([log
                   for log in wide_logs
                   if occursin("Common ThorAxe PID species list is small",
        String(log.message))])
end
