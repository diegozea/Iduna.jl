@testset "cached ThorAxe MSA result reuse" begin
    mktempdir() do tmp
        source = joinpath(tmp, "cached_input")
        ensembl = joinpath(source, "Ensembl")
        mkpath(ensembl)
        for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
            write(joinpath(ensembl, file), "$(file)\n")
        end

        workdir = joinpath(tmp, "work")
        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1")
        input_dir = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
            specieslist = "homo_sapiens",
            cached_input_dir = source,
            orthology = "1:1",
            overwrite = false)
        metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
            input_dir, target, [10.0];
            sample_count = 0,
            sample_fraction = 0.8,
            sample_seed = UInt64(7),
            requested_sample_seed = 7,
            sampling_strategy = :common,
            effective_specieslist = "homo_sapiens",
            orthology = "1:1",
            specieslist_filter = true,
            biomart_datasets_filter = true)

        paths = Iduna.ThorAxeMSA._pid_sample_paths(workdir, 10.0, 0)
        write(joinpath(tmp, "seed.fasta"), ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")
        seed_msa = Iduna.ThorAxeMSA.read_file(
            joinpath(tmp, "seed.fasta"), Iduna.ThorAxeMSA.FASTA)
        Iduna.Utils.set_s_exon_annotations!(seed_msa, "00", ['0' => "1_0"])
        mkpath(dirname(paths.fasta_path))
        Iduna.ThorAxeMSA.write_file(paths.fasta_path, seed_msa, Iduna.ThorAxeMSA.FASTA)
        Iduna.ThorAxeMSA.write_file(
            paths.stockholm_path, seed_msa, Iduna.ThorAxeMSA.Stockholm)
        Iduna.ThorAxeMSA._write_candidate_sample_inputs(
            paths, seed_msa, ["homo_sapiens", "mus_musculus"], [1, 2])
        thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(workdir, 10.0, 0)
        mkpath(thoraxe_dir)
        write(joinpath(thoraxe_dir, "path_table.csv"),
            "TranscriptIDCluster,Path\nENST00000000001,start/1_0/stop\n")

        candidate = (;
            msa = seed_msa,
            fasta_path = paths.fasta_path,
            stockholm_path = paths.stockholm_path,
            sequence_fasta = paths.sequence_fasta,
            species_file = paths.species_file,
            workdir)
        row = Iduna.ThorAxeMSA._candidate_summary_row(
            target, candidate, 10.0, 1,
            (; eligible = true, status = "ok", issue = missing);
            sample_count = metadata.pid_sample_count,
            sample_fraction = metadata.pid_sample_fraction,
            sample_seed = metadata.pid_sample_seed,
            metadata)
        summary_path = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir),
            "candidate_summary.csv")
        Iduna.ThorAxeMSA._summarize_candidate_scores([row], summary_path)
        summary_df = Iduna.ThorAxeMSA._candidate_summary_dataframe(summary_path)
        @test only(summary_df.stockholm_path) ==
              relpath(paths.stockholm_path, workdir)
        @test only(summary_df.fasta_path) == relpath(paths.fasta_path, workdir)
        @test only(summary_df.sequence_fasta) ==
              relpath(paths.sequence_fasta, workdir)
        @test only(summary_df.species_file) == relpath(paths.species_file, workdir)
        @test only(summary_df.scores_path) ==
              relpath(Iduna.ThorAxeMSA._pid_scores_path(workdir, 10.0), workdir)
        seed = Iduna.SeedSelection(;
            pid = 10.0,
            epli = missing,
            stockholm_path = paths.stockholm_path,
            fasta_path = paths.fasta_path,
            summary_path)
        Iduna.ThorAxeMSA._mark_selected_candidates!(summary_path, [seed])
        @test Iduna.ThorAxeMSA._cached_selected_seeds(
            summary_path, workdir, metadata) !== nothing

        stale_seed = Iduna.ThorAxeMSA.read_file(
            joinpath(tmp, "seed.fasta"), Iduna.ThorAxeMSA.FASTA)
        Iduna.ThorAxeMSA.write_file(
            paths.stockholm_path, stale_seed, Iduna.ThorAxeMSA.Stockholm)
        @test Iduna.ThorAxeMSA._cached_selected_seeds(
            summary_path, workdir, metadata) === nothing
        Iduna.ThorAxeMSA.write_file(
            paths.stockholm_path, seed_msa, Iduna.ThorAxeMSA.Stockholm)

        reuse_logs,
        result = Test.collect_test_logs() do
            Iduna.ThorAxeMSA.build_thoraxe_msa(
                target, workdir;
                pid_thresholds = [10.0],
                specieslist = "homo_sapiens",
                cached_thoraxe_input_dir = source,
                pid_sample_count = 0,
                pid_sample_seed = 7)
        end
        @test any(log -> log.message == "Reusing cached ThorAxe MSA candidates.",
            reuse_logs)
        @test !any(
            log -> log.message in ("Resolving ThorAxe species filters.",
                "Preparing ThorAxe transcript_query input."),
            reuse_logs)
        resolved_species_log = only([log
                                     for log in reuse_logs
                                     if log.message ==
                                        "Resolved ThorAxe transcript_query species."])
        resolved_species_kwargs = Dict(resolved_species_log.kwargs)
        @test resolved_species_kwargs[:n_requested_species] == 1
        @test resolved_species_kwargs[:n_effective_species] == 1
        @test !haskey(resolved_species_kwargs, :requested_specieslist_preview)
        @test !haskey(resolved_species_kwargs, :effective_specieslist_preview)
        @test result.input_dir == input_dir
        @test result.pid_sample_seed == UInt64(7)
        @test result.pid_sample_count == 0
        @test result.status === :ok
        @test isempty(result.warnings)
        @test only(result.seeds).pid == 10.0
        @test result.baseline_fastas == [paths.fasta_path]
        @test result.baseline_stockholms == [paths.stockholm_path]
        @test result.sequence_fastas == [paths.sequence_fasta]
        @test result.species_files == [paths.species_file]
        @test result.thoraxe_dirs == [thoraxe_dir]

        function write_cached_test_thoraxe_dir(thoraxe_root::AbstractString)
            msa_dir = joinpath(thoraxe_root, "msa")
            phylosofs_dir = joinpath(thoraxe_root, "phylosofs")
            mkpath(msa_dir)
            mkpath(phylosofs_dir)
            write(joinpath(thoraxe_root, "path_table.csv"),
                "TranscriptIDCluster,Path\nENST00000000001,start/1_0/stop\n")
            write(joinpath(thoraxe_root, "s_exon_table.csv"),
                "GeneID,Species,TranscriptIDCluster,S_exonID,S_exon_Sequence\n" *
                "ENSG00000000001,homo_sapiens,ENST00000000001,1_0,AA\n" *
                "ORTHO1,mus_musculus,ORTHO1,1_0,AB\n")
            write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
                ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")
            write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\ta\n")
            write(joinpath(phylosofs_dir, "transcripts.pir"),
                ">P1;ENSG00000000001 ENST00000000001 a\naa\nAA*\n")
            return thoraxe_root
        end

        scored_workdir = joinpath(tmp, "scored_work")
        scored_input = Iduna.ThorAxeMSA._ensure_transcript_query(
            target, scored_workdir;
            cached_input_dir = source,
            overwrite = false)
        scored_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
            scored_input, target, [20.0, 30.0];
            sample_count = 0,
            sample_fraction = 0.8,
            sample_seed = UInt64(9),
            requested_sample_seed = 9,
            sampling_strategy = :common,
            effective_specieslist = nothing,
            orthology = "1:1",
            specieslist_filter = false,
            biomart_datasets_filter = false)
        scored_summary = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(scored_workdir),
            "candidate_summary.csv")
        mkpath(dirname(scored_summary))
        scored_candidate = (;
            msa = seed_msa,
            fasta_path = "unselected.fasta",
            stockholm_path = "unselected.sto",
            sequence_fasta = "unselected_sequences.fasta",
            species_file = "unselected_species.txt",
            workdir = scored_workdir)
        scored_rows = [
            Iduna.ThorAxeMSA._candidate_summary_row(
                target, scored_candidate, 20.0, 1,
                (; eligible = true, status = "ok", issue = missing);
                sample_count = scored_metadata.pid_sample_count,
                sample_fraction = scored_metadata.pid_sample_fraction,
                sample_seed = scored_metadata.pid_sample_seed,
                metadata = scored_metadata),
            Iduna.ThorAxeMSA._candidate_summary_row(
                target, scored_candidate, 30.0, 2,
                (; eligible = true, status = "ok", issue = missing);
                sample_count = scored_metadata.pid_sample_count,
                sample_fraction = scored_metadata.pid_sample_fraction,
                sample_seed = scored_metadata.pid_sample_seed,
                metadata = scored_metadata)
        ]
        Iduna.ThorAxeMSA._summarize_candidate_scores(scored_rows, scored_summary)
        for pid in (20.0, 30.0)
            write_cached_test_thoraxe_dir(
                Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(scored_workdir, pid, 0))
        end
        scored_logs,
        scored_result = Test.collect_test_logs() do
            Iduna.ThorAxeMSA.build_thoraxe_msa(target, scored_workdir;
                pid_thresholds = [20.0, 30.0],
                specieslist = "all",
                cached_thoraxe_input_dir = source,
                specieslist_filter = false,
                biomart_datasets_filter = false,
                pid_sample_count = 0,
                pid_sample_seed = 9)
        end
        build_log = only([log
                          for log in scored_logs
                          if log.message == "Building ThorAxe MSA candidates."])
        @test !haskey(Dict(build_log.kwargs), :pid_sample_fraction)
        @test scored_result.status === :ok
        @test scored_result.pid_sample_count == 0
        @test [seed.pid for seed in scored_result.seeds] == [20.0, 30.0]
        @test length(scored_result.baseline_stockholms) == 2
        @test all(isfile, scored_result.baseline_stockholms)
        @test isfile(scored_result.pid_summary)
    end
end
