using Iduna, Test

@testset "ThorAxe symlink virtual environment provenance" begin
    msa = Iduna.ThorAxeMSA
    conda = msa.CondaPkg
    installed = conda.withenv() do
        python = Sys.which("python")
        prefix = strip(read(`$python -B -c "import sys; print(sys.prefix)"`, String))
        source = strip(read(
            `$python -B -c "import pathlib, thoraxe; print(pathlib.Path(thoraxe.__file__).parent)"`,
            String))
        scripts = Dict(name => read(Sys.which(name), String)
        for name in ("thoraxe", "transcript_query"))
        (; python, prefix, source, scripts)
    end

    mktempdir() do root
        venv = joinpath(realpath(root), "venv")
        run(`$(installed.python) -m venv --without-pip --system-site-packages --symlinks $venv`)
        python = joinpath(venv, "bin", "python")
        @test islink(python)
        @test realpath(python) == realpath(installed.python)
        site = strip(read(`$python -B -c "import sysconfig; print(sysconfig.get_path('purelib'))"`, String))
        package = joinpath(site, "thoraxe")
        cp(installed.source, package)
        rm(joinpath(package, "__pycache__"); recursive = true, force = true)
        function write_cli(name, interpreter = python)
            lines = split(installed.scripts[name], '\n')
            lines[1] = "#!" * interpreter
            path = joinpath(venv, "bin", name)
            write(path, join(lines, '\n'))
            chmod(path, 0o755)
            return path
        end
        cli = write_cli("thoraxe")
        write_cli("transcript_query")

        msa._with_thoraxe_lock() do
            original_backend = conda.STATE.backend
            try
                conda.STATE.backend = :Null
                withenv("PATH" => dirname(python) * ":" * ENV["PATH"],
                    "PYTHONPATH" => nothing, "PYTHONDONTWRITEBYTECODE" => "1") do
                    identities = Any[]
                    for version in ("0.8.3+venv_first", "0.8.3+venv_second_revision")
                        write(joinpath(package, "version.py"), "__version__ = \"$version\"\n")
                        observed = strip(read(`$cli --version`, String))
                        identity = msa._transcript_query_tool_identity()
                        @test occursin(version, observed)
                        @test identity !== nothing
                        identity === nothing && continue
                        @test identity.metadata["thoraxe_version"] == version
                        @test identity.metadata["python_executable"] == python
                        @test get(identity.metadata, "python_prefix", nothing) == venv
                        @test get(identity.metadata, "python_base_prefix", nothing) ==
                              installed.prefix
                        @test site in identity.metadata["import_paths"]
                        push!(identities, identity)
                    end
                    @test length(identities) == 2
                    if length(identities) == 2
                        before, after = identities
                        @test before.python == after.python
                        @test before.metadata["source_sha256"] !=
                              after.metadata["source_sha256"]
                        @test before != after

                        # Real tool identities drive the production cache decision;
                        # the bundle producer is local to avoid Ensembl requests.
                        target = Iduna.ResolvedTarget(;
                            input_id = "ENST", input_kind = :ensembl_transcript,
                            ensembl_gene_id = "ENSG", transcript_id = "ENST")
                        calls = Ref(0)
                        runner = command -> begin
                            calls[] += 1
                            ensembl = mkpath(joinpath(pwd(), last(command.exec), "Ensembl"))
                            for file in msa._REQUIRED_ENSEMBL_FILES
                                write(joinpath(ensembl, file), "fixture\n")
                            end
                            nothing
                        end
                        workdir = joinpath(root, "work")
                        ensure(identity) = msa._ensure_transcript_query(target, workdir;
                            max_retries = 1, transcript_query_runner = runner,
                            _implementation = (; source = "venv-regression"),
                            _tool_identity = identity)
                        ensure(before)
                        ensure(before)
                        @test calls[] == 1
                        ensure(after)
                        @test calls[] == 2
                        @test Iduna.Utils._read_stage_state(
                            msa._thoraxe_input_stage_dir(workdir))["action"] == "rebuild"
                    end

                    other_venv = joinpath(root, "other_venv")
                    run(`$(installed.python) -m venv --without-pip --system-site-packages --symlinks $other_venv`)
                    other_python = joinpath(other_venv, "bin", "python")
                    @test realpath(other_python) == realpath(python)
                    write_cli("transcript_query", other_python)
                    @test_logs (:warn, r"ThorAxe tool identity unavailable") begin
                        @test msa._transcript_query_tool_identity() === nothing
                    end
                end
            finally
                conda.STATE.backend = original_backend
            end
        end
    end
end
