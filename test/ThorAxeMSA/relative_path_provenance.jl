using Iduna, Test

@testset "ThorAxe executable lookup across working directories" begin
    msa = Iduna.ThorAxeMSA
    conda = msa.CondaPkg
    python = conda.withenv() do
        Sys.which("python")
    end
    separator = Sys.iswindows() ? ';' : ':'
    target = Iduna.ResolvedTarget(;
        input_id = "ENST", input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG", transcript_id = "ENST")

    mktempdir() do tmp
        root = realpath(tmp)
        caller = mkpath(joinpath(root, "caller"))
        discovered_bin = mkpath(joinpath(root, "bin"))
        workdir = mkpath(joinpath(root, "work"))
        actual_bin = mkpath(joinpath(workdir, "bin"))
        executions = joinpath(root, "executions.txt")
        function write_producer(path, label)
            body = "#!" * python * "\nimport pathlib, sys\n"
            body *= "dest = pathlib.Path(sys.argv[-1]) / 'Ensembl'\n"
            body *= "dest.mkdir(parents=True, exist_ok=True)\n"
            for file in msa._REQUIRED_ENSEMBL_FILES
                body *= "(dest / " * repr(file) * ").write_text(" * repr(label * "\n") *
                        ")\n"
            end
            body *= "(dest / 'executed-path.txt').write_text(str(pathlib.Path(sys.argv[0]).resolve()))\n"
            body *= "with pathlib.Path(" * repr(executions) * ").open('a') as log:\n"
            body *= "    log.write(" * repr(label * "\n") * ")\n"
            write(path, body)
            chmod(path, 0o755)
            return path
        end
        for name in ("thoraxe", "transcript_query")
            write_producer(joinpath(discovered_bin, name), "caller-A")
            write_producer(joinpath(actual_bin, name), "actual-B")
        end
        call_count() = isfile(executions) ? countlines(executions) : 0
        content(input) = read(joinpath(input, "Ensembl", first(msa._REQUIRED_ENSEMBL_FILES)), String)
        action(work) = Iduna.Utils._read_stage_state(msa._thoraxe_input_stage_dir(work))["action"]
        ensure(work; overwrite = false) = msa._ensure_transcript_query(target, work;
            max_retries = 1, overwrite)
        absolute_path = join((discovered_bin, dirname(python)), separator)
        relative_path = join((joinpath("..", "bin"), dirname(python)), separator)

        msa._with_thoraxe_lock() do
            original_backend = conda.STATE.backend
            original_environment = copy(ENV)
            try
                conda.STATE.backend = :Null
                cd(caller) do
                    withenv("PATH" => relative_path, "PYTHONPATH" => nothing,
                        "PYTHONDONTWRITEBYTECODE" => "1") do
                        @test Sys.which("transcript_query") ==
                              joinpath(discovered_bin, "transcript_query")
                        @test_logs (:warn, r"ThorAxe tool identity unavailable") begin
                            @test msa._transcript_query_tool_identity() === nothing
                        end
                        # The default subprocess executes B after production changes cwd.
                        input = ensure(workdir)
                        @test content(input) == "actual-B\n"
                        @test read(joinpath(input, "Ensembl", "executed-path.txt"), String) ==
                              joinpath(actual_bin, "transcript_query")
                        @test call_count() == 1
                        write_producer(joinpath(actual_bin, "transcript_query"), "actual-C")
                        ensure(workdir)
                        @test content(input) == "actual-C\n"
                        @test action(workdir) == "rebuild"
                        @test call_count() == 2
                        ensure(workdir)
                        @test content(input) == "actual-C\n"
                        @test action(workdir) == "rebuild"
                        @test call_count() == 3
                        cached_content = content(input)
                        ensure(workdir; overwrite = true)
                        @test content(input) == cached_content
                        @test call_count() == 4
                    end

                    withenv("PATH" => absolute_path, "PYTHONPATH" => nothing,
                        "PYTHONDONTWRITEBYTECODE" => "1") do
                        before = msa._transcript_query_tool_identity()
                        @test before !== nothing
                        absolute_work = joinpath(root, "absolute_work")
                        input = ensure(absolute_work)
                        @test content(input) == "caller-A\n"
                        @test call_count() == 5
                        ensure(absolute_work)
                        @test action(absolute_work) == "reuse"
                        @test call_count() == 5
                        write_producer(joinpath(discovered_bin, "transcript_query"), "caller-D")
                        @test msa._transcript_query_tool_identity() != before
                        ensure(absolute_work)
                        @test content(input) == "caller-D\n"
                        @test action(absolute_work) == "rebuild"
                        @test call_count() == 6
                    end

                    # Empty entries and currently absent relative directories can also
                    # select different commands after a later working-directory change.
                    for path in (".", "", nothing,
                        "missing" * separator * absolute_path,
                        separator * absolute_path, absolute_path * separator,
                        absolute_path * separator * separator * dirname(python))
                        withenv("PATH" => path) do
                            @test_logs (:warn, r"ThorAxe tool identity unavailable") begin
                                @test msa._transcript_query_tool_identity() === nothing
                            end
                        end
                    end
                end
            finally
                conda.STATE.backend = original_backend
            end
            @test copy(ENV) == original_environment
            @test conda.STATE.backend == original_backend
        end
    end
end
