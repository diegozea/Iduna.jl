@testset "cached ThorAxe MSA result reuse" begin
    mktempdir() do tmp
        provenance = (; implementation = (; source = "fixture-v1"),
            python = (; source = "python-fixture-v1"),
            aligner = (; sha256 = "aligner-fixture-v1"), hhsuite = nothing)
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
            overwrite = false,
            _implementation = provenance.implementation)
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
            biomart_datasets_filter = true,
            provenance)

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

        stage_identity = Iduna.ThorAxeMSA._thoraxe_msa_identity(metadata, workdir)
        Iduna.ThorAxeMSA._write_thoraxe_msa_state(workdir, summary_path, [seed],
            :done, stage_identity; action = :run)
        reuse_logs,
        result = Test.collect_test_logs() do
            Iduna.ThorAxeMSA.build_thoraxe_msa(
                target, workdir;
                pid_thresholds = [10.0],
                specieslist = "homo_sapiens",
                cached_thoraxe_input_dir = source,
                pid_sample_count = 0,
                pid_sample_seed = 7,
                _provenance = provenance)
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

        @testset "public ThorAxe relative workdirs and relocation" begin
            for path_kind in (:absolute, :relative, :dot)
                @testset "$path_kind" begin
                    relocated = joinpath(tmp, "moved_$(path_kind)")
                    cp(workdir, relocated)
                    cd(path_kind === :dot ? relocated : tmp) do
                        requested = path_kind === :relative ? basename(relocated) :
                                    path_kind === :dot ? "." : relocated
                        logs, moved = Test.collect_test_logs() do
                            Iduna.ThorAxeMSA.build_thoraxe_msa(target, requested;
                                pid_thresholds = [10.0], specieslist = "homo_sapiens",
                                cached_thoraxe_input_dir = source, pid_sample_count = 0,
                                pid_sample_seed = 7, _provenance = provenance)
                        end
                        @test any(
                            log -> log.message ==
                                   "Reusing cached ThorAxe MSA candidates.", logs)
                        @test isabspath(moved.input_dir)
                        @test realpath(moved.input_dir) ==
                              realpath(joinpath(relocated, "thoraxe_input"))
                        @test all(isfile, moved.baseline_stockholms)
                        @test all(isfile, moved.sequence_fastas)
                        @test all(isfile, moved.species_files)
                        state = Iduna.Utils._read_stage_state(
                            Iduna.ThorAxeMSA._thoraxe_msa_stage_dir(relocated))
                        @test state["action"] == "reuse"
                        @test all(!isabspath(path)
                        for path in keys(state["output_digests"]))
                    end
                end
            end
        end

        # Independently edit each required artifact after a sound completion.
        # Metadata-only and candidate-file fallbacks must not turn a miss into a hit.
        for file in (paths.stockholm_path, paths.fasta_path, paths.sequence_fasta,
            paths.species_file, joinpath(thoraxe_dir, "path_table.csv"), summary_path)
            original = read(file, String)
            write(file, replace(original, first(original) => '!'; count = 1))
            cache = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(workdir,
                summary_path, metadata, stage_identity; overwrite = false)
            @test !cache.cache.reusable
            @test Iduna.ThorAxeMSA._maybe_cached_thoraxe_msa(input_dir, summary_path,
                target, workdir, metadata, stage_identity, cache, nothing, nothing;
                overwrite = false, pid_sample_count = 0, pid_sample_fraction = 0.8) ===
                  nothing
            write(file, original)
        end
        stage_path = Iduna.Utils._stage_state_path(Iduna.ThorAxeMSA._thoraxe_msa_stage_dir(workdir))
        valid_manifest = read(stage_path, String)
        write(stage_path, "{}")
        stale = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(workdir,
            summary_path, metadata, stage_identity; overwrite = false)
        @test !stale.cache.reusable
        @test Iduna.ThorAxeMSA._maybe_cached_thoraxe_msa(input_dir, summary_path,
            target, workdir, metadata, stage_identity, stale, nothing, nothing;
            overwrite = false, pid_sample_count = 0, pid_sample_fraction = 0.8) === nothing
        write(stage_path, valid_manifest)
        incomplete_manifest = Iduna.Utils.JSON.parse(valid_manifest)
        delete!(incomplete_manifest["required_outputs"], "seed_fastas")
        delete!(incomplete_manifest["output_digests"], relpath(paths.fasta_path, workdir))
        Iduna.Utils.write_json(stage_path, incomplete_manifest)
        incomplete = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(workdir,
            summary_path, metadata, stage_identity; overwrite = false)
        @test incomplete.cache.reusable
        @test Iduna.ThorAxeMSA._maybe_cached_thoraxe_msa(input_dir, summary_path,
            target, workdir, metadata, stage_identity, incomplete, nothing, nothing;
            overwrite = false, pid_sample_count = 0, pid_sample_fraction = 0.8) === nothing
        write(stage_path, valid_manifest)

        rm(paths.s_exon_blocks_tsv)
        cached_again = Iduna.ThorAxeMSA._cached_selected_seeds(summary_path, workdir, metadata)
        @test cached_again !== nothing
        @test isfile(paths.s_exon_blocks_tsv)

        auto_metadata = merge(metadata, (; requested_pid_sample_seed = nothing))
        auto_identity = Iduna.ThorAxeMSA._thoraxe_msa_identity(auto_metadata, workdir)
        Iduna.ThorAxeMSA._write_thoraxe_msa_state(workdir, summary_path, [seed],
            :done, auto_identity; action = :run)
        auto_reuse = Iduna.ThorAxeMSA.build_thoraxe_msa(target, workdir;
            pid_thresholds = [10.0], specieslist = "homo_sapiens",
            cached_thoraxe_input_dir = source, pid_sample_count = 0,
            _provenance = provenance)
        @test auto_reuse.pid_sample_seed == UInt64(7)

        original_sto = read(paths.stockholm_path, String)
        write(paths.stockholm_path, "not Stockholm")
        @test !Iduna.ThorAxeMSA._selected_artifacts_exist(seed, workdir)
        @test_throws ErrorException Iduna.ThorAxeMSA._write_thoraxe_msa_state(workdir,
            summary_path, [seed], :done, auto_identity; action = :run)
        write(paths.stockholm_path, original_sto)
        original_fasta = read(paths.fasta_path, String)
        write(paths.fasta_path, ">ENSG00000000001\n")
        @test !Iduna.ThorAxeMSA._selected_artifacts_exist(seed, workdir)
        @test_throws ErrorException Iduna.ThorAxeMSA._write_thoraxe_msa_state(workdir,
            summary_path, [seed], :done, auto_identity; action = :run)
        write(paths.fasta_path, original_fasta)
        @test_throws ErrorException Iduna.ThorAxeMSA._write_thoraxe_msa_state(workdir,
            summary_path, [seed], :done, merge(auto_identity, (; pid_path_contract = 1));
            action = :run)
        alternative_fasta = joinpath(workdir, "other_seed.fasta")
        cp(paths.fasta_path, alternative_fasta)
        different_path = Iduna.SeedSelection(; pid = seed.pid, epli = seed.epli,
            stockholm_path = seed.stockholm_path, fasta_path = alternative_fasta,
            summary_path)
        @test_throws ErrorException Iduna.ThorAxeMSA._write_thoraxe_msa_state(workdir,
            summary_path, [different_path], :done, auto_identity; action = :run)

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
            overwrite = false,
            _implementation = provenance.implementation)
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
            biomart_datasets_filter = false,
            provenance)
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
            identity = Iduna.ThorAxeMSA._thoraxe_msa_identity(scored_metadata, scored_workdir)
            cache = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(scored_workdir,
                scored_summary, scored_metadata, identity; overwrite = false)
            @test !cache.cache.reusable
            prepared = Iduna.ThorAxeMSA._prepare_thoraxe_msa_stage!(scored_workdir,
                scored_summary, identity, cache; overwrite = false)
            @test prepared.force_pid_rerun
            fake_thoraxe = (
                input, run_root; kwargs...) -> write_cached_test_thoraxe_dir(joinpath(run_root, "thoraxe"))
            filters = (; species_filter = (; warnings = String[]),
                biomart_filter = (; warnings = String[]))
            Iduna.ThorAxeMSA._run_thoraxe_msa_stage!(target, scored_input,
                scored_workdir, scored_summary, [20.0, 30.0], nothing,
                scored_metadata, identity, filters, prepared;
                pid_sample_count = 0,
                pid_sample_fraction = 0.8, sample_seed = UInt64(9),
                sampling_strategy = :common, thoraxe_fn = fake_thoraxe)
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
        rerun = Iduna.ThorAxeMSA.build_thoraxe_msa(target, scored_workdir;
            pid_thresholds = [20.0, 30.0], specieslist = "all",
            cached_thoraxe_input_dir = source, specieslist_filter = false,
            biomart_datasets_filter = false, pid_sample_count = 0, pid_sample_seed = 9,
            _provenance = provenance)
        @test [seed.pid for seed in rerun.seeds] == [20.0, 30.0]

        # A rejected stage must attempt a rebuild and retain failed state if its
        # output directory is blocked. This fails locally before invoking ThorAxe.
        blocked_workdir = joinpath(tmp, "blocked_work")
        mkpath(blocked_workdir)
        blocked_output = Iduna.ThorAxeMSA._thoraxe_msa_dir(blocked_workdir)
        write(blocked_output, "blocking file")
        @test_throws Base.IOError Iduna.ThorAxeMSA.build_thoraxe_msa(
            target, blocked_workdir;
            pid_thresholds = [20.0], specieslist = "all",
            cached_thoraxe_input_dir = source, specieslist_filter = false,
            biomart_datasets_filter = false, pid_sample_count = 0, pid_sample_seed = 9,
            _provenance = provenance)
        failed = Iduna.Utils._read_stage_state(
            Iduna.ThorAxeMSA._thoraxe_msa_stage_dir(blocked_workdir))
        @test failed["status"] == "failed"
        @test failed["output_digests"] === nothing
        @test read(blocked_output, String) == "blocking file"

        rebuilds = Ref(0)
        fixture_runner = (args...; kwargs...) -> begin
            rebuilds[] += 1
            fake = (input, root; kwargs...) -> write_cached_test_thoraxe_dir(joinpath(root, "thoraxe"))
            Iduna.ThorAxeMSA._run_thoraxe_msa_stage!(args...; kwargs..., thoraxe_fn = fake)
        end
        function build_fixture(path, facts; query = target, species = "all")
            return Iduna.ThorAxeMSA.build_thoraxe_msa(query, path;
                pid_thresholds = [20.0, 30.0], specieslist = species,
                cached_thoraxe_input_dir = source, specieslist_filter = false,
                biomart_datasets_filter = false, pid_sample_count = 0, pid_sample_seed = 9,
                _provenance = facts, _stage_runner = fixture_runner)
        end
        # Every producer change reaches the real public rebuild path. The fixture
        # replaces only the external biological tool, retaining assembly and writes.
        for field in (:implementation, :python, :aligner)
            changed = merge(provenance, NamedTuple{(field,)}(((; revision = "changed"),)))
            path = joinpath(tmp, "changed_$(field)")
            cp(scored_workdir, path)
            before = rebuilds[]
            build_fixture(path, changed)
            @test rebuilds[] == before + 1
            build_fixture(path, changed)
            @test rebuilds[] == before + 1

            unknown = merge(changed, NamedTuple{(field,)}((nothing,)))
            build_fixture(path, unknown)
            build_fixture(path, unknown)
            @test rebuilds[] == before + 3
        end

        reference = joinpath(tmp, "reference.fasta")
        write(reference, ">P1\nAA\n")
        reference_target = Iduna.ResolvedTarget(;
            input_id = target.input_id, input_kind = target.input_kind,
            ensembl_gene_id = target.ensembl_gene_id, transcript_id = target.transcript_id,
            uniprot_id = "P1", uniprot_sequence_path = reference)
        reference_work = joinpath(tmp, "reference_work")
        cp(scored_workdir, reference_work)
        before = rebuilds[]
        build_fixture(reference_work, provenance; query = reference_target)
        build_fixture(reference_work, provenance; query = reference_target)
        @test rebuilds[] == before + 1
        write(reference, ">P1\nAC\n")
        build_fixture(reference_work, provenance; query = reference_target)
        @test rebuilds[] == before + 2

        species_file = joinpath(tmp, "requested_species.txt")
        write(species_file, "homo_sapiens\n")
        species_work = joinpath(tmp, "species_work")
        cp(scored_workdir, species_work)
        before = rebuilds[]
        build_fixture(species_work, provenance; species = species_file)
        build_fixture(species_work, provenance; species = species_file)
        @test rebuilds[] == before + 1
        write(species_file, "mus_musculus\n")
        build_fixture(species_work, provenance; species = species_file)
        @test rebuilds[] == before + 2
    end
end
