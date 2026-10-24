@testset "ThorAxe input producer provenance" begin
    mktempdir() do tmp
        target = Iduna.ResolvedTarget(;
            input_id = "ENST", input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG", transcript_id = "ENST")
        implementation = (; source_sha256 = "test-source-a")
        python = (; executable_sha256 = "test-python-a", source_sha256 = "test-thoraxe-a")
        calls = Ref(0)
        runner = command -> begin
            calls[] += 1
            write_test_ensembl_bundle(joinpath(pwd(), last(command.exec)))
            nothing
        end
        workdir = joinpath(tmp, "queried")
        function ensure(; implementation = implementation, python = python,
                specieslist = nothing)
            return Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
                specieslist, max_retries = 1, transcript_query_runner = runner,
                _implementation = implementation, _tool_identity = python)
        end

        input_dir = ensure()
        @test calls[] == 1
        @test ensure() == input_dir
        @test calls[] == 1
        @test Iduna.Utils._read_stage_state(
            Iduna.ThorAxeMSA._thoraxe_input_stage_dir(workdir))["action"] == "reuse"
        changed_python = merge(python, (; source_sha256 = "test-thoraxe-b"))
        ensure(; python = changed_python)
        @test calls[] == 2
        ensure(; implementation = (; source_sha256 = "test-source-b"),
            python = changed_python)
        @test calls[] == 3

        # Repeating the same unknown identity must not create a reusable key.
        ensure(; python = nothing)
        ensure(; python = nothing)
        @test calls[] == 5
        ensure(; implementation = nothing)
        ensure(; implementation = nothing)
        @test calls[] == 7

        species_file = joinpath(tmp, "species.txt")
        write(species_file, "homo_sapiens\n")
        ensure(; specieslist = species_file)
        @test calls[] == 8
        ensure(; specieslist = species_file)
        @test calls[] == 8
        write(species_file, "mus_musculus\n")
        ensure(; specieslist = species_file)
        @test calls[] == 9

        source = write_test_ensembl_bundle(joinpath(tmp, "copied_source"))
        copied_workdir = joinpath(tmp, "copied")
        copy_input() = Iduna.ThorAxeMSA._ensure_transcript_query(target, copied_workdir;
            cached_input_dir = source, _implementation = implementation,
            _tool_identity = nothing,
            transcript_query_runner = command -> error("Copied input must not run Python."))
        copied_input = copy_input()
        @test copy_input() == copied_input
        copied_state = Iduna.Utils._read_stage_state(
            Iduna.ThorAxeMSA._thoraxe_input_stage_dir(copied_workdir))
        @test copied_state["action"] == "reuse"
        source_fasta = joinpath(source, "Ensembl", "sequences.fasta")
        write(source_fasta, ">changed\nACDE\n")
        copy_input()
        @test read(joinpath(copied_input, "Ensembl", "sequences.fasta"), String) ==
              ">changed\nACDE\n"
        @test Iduna.Utils._read_stage_state(
            Iduna.ThorAxeMSA._thoraxe_input_stage_dir(copied_workdir))["action"] ==
              "rebuild"
    end
end

@testset "ThorAxe tool identity failure boundaries" begin
    thoraxe_msa = Iduna.ThorAxeMSA
    @test thoraxe_msa._inspect_thoraxe_identity(() -> 7) == 7
    @test_throws InterruptException thoraxe_msa._inspect_thoraxe_identity(
        () -> throw(InterruptException()))
    failed = @test_logs (:warn, r"tool identity unavailable") thoraxe_msa._inspect_thoraxe_identity(
        () -> error("expected identity failure"))
    @test failed === nothing

    mktempdir() do tmp
        interpreter = joinpath(tmp, "python")
        write(interpreter,
            "#!/bin/sh\nprintf '%s\\n' " *
            "'{\"python_executable\":\"/bin/sh\",\"source_sha256\":\"test-source\"}'\n")
        chmod(interpreter, 0o755)
        function write_cli(path, executable = interpreter)
            write(path, "#!$(executable)\n")
            chmod(path, 0o755)
            return path
        end
        write_cli(joinpath(tmp, "thoraxe"))
        query = write_cli(joinpath(tmp, "transcript_query"))
        thoraxe = thoraxe_msa.ThorAxe
        conda = thoraxe.CondaPkg

        # Isolate discovery from installed Python; restore shared wrapper state even
        # if an assertion fails. The shell fixture returns probe JSON without Python.
        thoraxe_msa._with_thoraxe_lock() do
            original_backend = conda.STATE.backend
            original_aligner = thoraxe.ALIGNER_PATH[]
            original_env = copy(ENV)
            original_frozen = conda.STATE.frozen
            try
                conda.STATE.backend = :Null
                thoraxe.ALIGNER_PATH[] = joinpath(tmp, "missing_aligner")
                withenv("PATH" => tmp) do
                    python_only = thoraxe_msa._thoraxe_tool_identity(; aligner = false)
                    @test python_only !== nothing
                    @test python_only.python.metadata["source_sha256"] == "test-source"
                    @test python_only.aligner === nothing
                    full = thoraxe_msa._thoraxe_tool_identity()
                    @test full.python == python_only.python
                    @test full.aligner === nothing

                    rm(query)
                    @test thoraxe_msa._thoraxe_tool_identity(; aligner = false) === nothing
                    write_cli(query, "/usr/bin/env python")
                    @test thoraxe_msa._thoraxe_tool_identity(; aligner = false) === nothing
                    other_interpreter = joinpath(tmp, "other_python")
                    cp(interpreter, other_interpreter)
                    write_cli(query, other_interpreter)
                    @test thoraxe_msa._thoraxe_tool_identity(; aligner = false) === nothing

                    other_dir = mkpath(joinpath(tmp, "other_directory"))
                    other_query = write_cli(joinpath(other_dir, "transcript_query"))
                    rm(query)
                    symlink(other_query, query)
                    @test thoraxe_msa._thoraxe_tool_identity(; aligner = false) === nothing
                    rm(query)
                    write_cli(query)
                    write(interpreter,
                        "#!/bin/sh\nprintf '%s\\n' " *
                        "'{\"python_executable\":\"/missing-test-python\"}'\n")
                    @test thoraxe_msa._thoraxe_tool_identity(; aligner = false) === nothing
                    write(interpreter, "#!/bin/sh\nexit 3\n")
                    @test thoraxe_msa._thoraxe_tool_identity(; aligner = false) === nothing
                end
            finally
                conda.STATE.backend = original_backend
                thoraxe.ALIGNER_PATH[] = original_aligner
            end
            env_restored = copy(ENV) == original_env
            @test env_restored
            @test conda.STATE.frozen == original_frozen
        end
    end
end

@testset "ThorAxe content and generated-seed identities" begin
    mktempdir() do tmp
        input_dir = write_test_ensembl_bundle(joinpath(tmp, "input"))
        reference_root = mkpath(joinpath(tmp, "references"))
        stage_root = mkpath(joinpath(tmp, "stage"))
        provenance = (;
            implementation = (; source_sha256 = "source-a"),
            python = (; source_sha256 = "python-a"),
            aligner = (; sha256 = "aligner-a"),
            hhsuite = (; sha256 = "hhsuite-a"))
        target(; kwargs...) = Iduna.ResolvedTarget(;
            input_id = "ENST", input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG", transcript_id = "ENST", uniprot_id = "P1",
            kwargs...)
        function metadata(reference_target; specieslist = nothing, provenance = provenance,
                sample_count = 2, seed = UInt64(7), requested_seed = nothing)
            return Iduna.ThorAxeMSA._candidate_run_metadata(
                input_dir, reference_target, [35.0];
                sample_count, sample_fraction = 0.5, sample_seed = seed,
                requested_sample_seed = requested_seed, sampling_strategy = :common,
                effective_specieslist = specieslist, orthology = "1:1",
                specieslist_filter = false, biomart_datasets_filter = false,
                provenance, workdir = stage_root)
        end
        absent = metadata(target())
        relative_target = target(;
            uniprot_sequence_path = "reference.fasta", workdir = reference_root)
        missing_reference = metadata(relative_target)
        @test absent.reference_identity != missing_reference.reference_identity
        reference = joinpath(reference_root, "reference.fasta")
        write(reference, ">P1\nACDE\n")
        present = metadata(relative_target)
        @test present.reference_identity != missing_reference.reference_identity
        @test present.reference_identity == metadata(target(;
            uniprot_sequence_path = reference)).reference_identity
        # The stage directory contains a decoy: target.workdir determines the reference.
        write(joinpath(stage_root, "reference.fasta"), ">P1\nWWWW\n")
        @test metadata(relative_target).reference_identity == present.reference_identity
        fallback_reference = metadata(target(; uniprot_sequence_path = "reference.fasta"))
        @test fallback_reference.reference_identity != present.reference_identity
        write(reference, ">P1\nACDF\n")
        changed_reference = metadata(relative_target)
        @test changed_reference.reference_identity != present.reference_identity

        species_file = joinpath(tmp, "species.txt")
        write(species_file, "homo_sapiens\n")
        first_species = metadata(relative_target; specieslist = species_file)
        write(species_file, "mus_musculus\n")
        second_species = metadata(relative_target; specieslist = species_file)
        @test first_species.effective_specieslist == second_species.effective_specieslist
        @test first_species.specieslist_fingerprint !=
              second_species.specieslist_fingerprint
        @test present.provenance == provenance

        canonical = Iduna.ThorAxeMSA._canonical_thoraxe_msa_metadata(present)
        @test canonical.pid_sample_seed === nothing
        @test canonical.provenance == present.provenance
        @test canonical.reference_identity == present.reference_identity
        @test Iduna.Utils._identity_hash(canonical) == Iduna.Utils._identity_hash(
            Iduna.ThorAxeMSA._canonical_thoraxe_msa_metadata(
            merge(present, (; pid_sample_seed = UInt64(99)))))
        for replacement in (
            (; provenance = merge(provenance, (; python = (; source_sha256 = "python-b")))),
            (; reference_identity = changed_reference.reference_identity),
            (; specieslist_fingerprint = second_species.specieslist_fingerprint))
            changed = Iduna.ThorAxeMSA._canonical_thoraxe_msa_metadata(
                merge(present, replacement))
            @test Iduna.Utils._identity_hash(changed) !=
                  Iduna.Utils._identity_hash(canonical)
        end
        explicit_seed = metadata(relative_target; requested_seed = 7)
        @test Iduna.ThorAxeMSA._canonical_thoraxe_msa_metadata(explicit_seed) ==
              explicit_seed

        # A digest-valid completed manifest is insufficient without known producers.
        summary = joinpath(stage_root, "summary.csv")
        write(summary, "pid\n35.0\n")
        function stage_cache(facts; sample_count = 2)
            run_metadata = metadata(relative_target; provenance = facts, sample_count)
            identity = Iduna.ThorAxeMSA._thoraxe_msa_identity(run_metadata, stage_root)
            stage_dir = Iduna.ThorAxeMSA._thoraxe_msa_stage_dir(stage_root)
            Iduna.Utils._write_stage_state(stage_dir;
                stage = "thoraxe_msa", stage_key = "thoraxe_msa", status = :done,
                identity, outputs = (; pid_summary = summary), workdir = stage_root)
            @test Iduna.Utils._classify_stage_state(stage_dir, identity,
                (; pid_summary = summary); stage_label = "test MSA").reusable
            return Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(stage_root, summary,
                run_metadata, identity; overwrite = false)
        end
        @test stage_cache(provenance).cache.reusable
        @test !stage_cache(nothing).cache.reusable
        for field in (:implementation, :python, :aligner, :hhsuite)
            unknown = merge(provenance, NamedTuple{(field,)}((nothing,)))
            @test !stage_cache(unknown).cache.reusable
        end
        @test stage_cache(merge(provenance, (; hhsuite = nothing));
            sample_count = 0).cache.reusable
    end
end
@testset "isolated installed tool provenance" begin
    script = joinpath(@__DIR__, "tool_provenance_process.jl")
    project = dirname(Base.active_project())
    # Preserve Pkg.test's coverage scope instead of instrumenting every dependency.
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project $script`
    process = run(pipeline(ignorestatus(command); stdout, stderr))
    @test (process.exitcode, process.termsignal) == (0, 0)
end
