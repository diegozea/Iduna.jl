mktempdir() do tmp
    stdout_log = joinpath(tmp, "stdout.log")
    stderr_log = joinpath(tmp, "stderr.log")
    quiet_terminal_path = joinpath(tmp, "quiet_terminal_stdout.log")
    open(quiet_terminal_path, "w") do terminal_io
        redirect_stdout(terminal_io) do
            Iduna.ThorAxeMSA._run_logged_command(
                `$(Base.julia_cmd()) --startup-file=no -e "println(\"ok\")"`,
                stdout_log,
                stderr_log)
        end
    end
    @test read(stdout_log, String) == "ok\n"
    @test isempty(read(quiet_terminal_path, String))
    @test isempty(read(stderr_log, String))

    live_stdout_path = joinpath(tmp, "live_process_stdout.log")
    live_stderr_path = joinpath(tmp, "live_process_stderr.log")
    open(live_stdout_path, "w") do live_stdout_io
        open(live_stderr_path, "w") do live_stderr_io
            redirect_stdout(live_stdout_io) do
                redirect_stderr(live_stderr_io) do
                    Iduna.ThorAxeMSA._run_logged_command(
                        `sh -c "printf live"`,
                        joinpath(tmp, "tee_stdout.log"),
                        joinpath(tmp, "tee_stderr.log");
                        live_stdout = true)
                end
            end
        end
    end
    @test isempty(read(live_stdout_path, String))
    @test read(live_stderr_path, String) == "live"
    @test read(joinpath(tmp, "tee_stdout.log"), String) == "live"
    @test isempty(read(joinpath(tmp, "tee_stderr.log"), String))

    query_runner = Iduna.ThorAxeMSA._transcript_query_thoraxe_runner(
        joinpath(tmp, "query_stdout.log"),
        joinpath(tmp, "query_stderr.log"))
    query_stdout_path = joinpath(tmp, "query_process_stdout.log")
    query_stderr_path = joinpath(tmp, "query_process_stderr.log")
    open(query_stdout_path, "w") do query_stdout_io
        open(query_stderr_path, "w") do query_stderr_io
            redirect_stdout(query_stdout_io) do
                redirect_stderr(query_stderr_io) do
                    query_runner(`sh -c "printf query-live"`)
                end
            end
        end
    end
    @test isempty(read(query_stdout_path, String))
    @test read(query_stderr_path, String) == "query-live"
    @test read(joinpath(tmp, "query_stdout.log"), String) == "query-live"
    @test isempty(read(joinpath(tmp, "query_stderr.log"), String))

    query_spinner_path = joinpath(tmp, "query_spinner_stderr.log")
    open(query_spinner_path, "w") do spinner_io
        redirect_stderr(spinner_io) do
            Iduna.ThorAxeMSA._run_logged_command(
                `sh -c "sleep 0.1; printf query-spinner"`,
                joinpath(tmp, "query_spinner_stdout.log"),
                joinpath(tmp, "query_spinner_command_stderr.log");
                live_stdout = true,
                progress_desc = "Running ThorAxe transcript_query: ",
                progress_enabled = true)
        end
    end
    query_spinner_text = read(query_spinner_path, String)
    @test Iduna.ThorAxeMSA._TRANSCRIPT_QUERY_SPINNER_INTERVAL_SECONDS ≈ 1 / 3
    @test occursin("Running ThorAxe transcript_query", query_spinner_text)
    @test !occursin("query-spinner", query_spinner_text)
    @test read(joinpath(tmp, "query_spinner_stdout.log"), String) ==
          "query-spinner"
    @test isempty(read(joinpath(tmp, "query_spinner_command_stderr.log"), String))

    withenv("CI" => "true", "GITHUB_ACTIONS" => nothing) do
        @test !Iduna.Utils._terminal_progress_enabled(IOBuffer())
    end
    withenv("CI" => nothing, "GITHUB_ACTIONS" => nothing) do
        @test !Iduna.Utils._terminal_progress_enabled(IOBuffer())
    end
    @test !Iduna.Utils._io_is_tty(IOContext(IOBuffer()))

    hidden_progress_path = joinpath(tmp, "hidden_progress.log")
    open(hidden_progress_path, "w") do progress_io
        Iduna.ThorAxeMSA._run_logged_command(
            `sh -c "printf hidden"`,
            joinpath(tmp, "hidden_stdout.log"),
            joinpath(tmp, "hidden_stderr.log");
            progress_desc = "Hidden progress: ",
            progress_output = progress_io)
    end
    @test isempty(read(hidden_progress_path, String))
    @test read(joinpath(tmp, "hidden_stdout.log"), String) == "hidden"

    spinner_path = joinpath(tmp, "spinner.log")
    open(spinner_path, "w") do progress_io
        Iduna.ThorAxeMSA._run_logged_command(
            `sh -c "sleep 0.2; printf spinner"`,
            joinpath(tmp, "spinner_stdout.log"),
            joinpath(tmp, "spinner_stderr.log");
            progress_desc = "Spinner progress: ",
            progress_output = progress_io,
            progress_enabled = true)
    end
    @test occursin("Spinner progress", read(spinner_path, String))
    @test read(joinpath(tmp, "spinner_stdout.log"), String) == "spinner"

    failed_spinner_path = joinpath(tmp, "failed_spinner.log")
    open(failed_spinner_path, "w") do progress_io
        @test_throws ErrorException Iduna.ThorAxeMSA._run_logged_command(
            `sh -c "sleep 0.1; exit 2"`,
            joinpath(tmp, "failed_spinner_stdout.log"),
            joinpath(tmp, "failed_spinner_stderr.log");
            progress_desc = "Failing progress: ",
            progress_output = progress_io,
            progress_enabled = true)
    end
    @test occursin("External command stopped",
        read(failed_spinner_path, String))

    held_stdout_log = joinpath(tmp, "held_stdout.log")
    held_stderr_log = joinpath(tmp, "held_stderr.log")
    started = time()
    Iduna.ThorAxeMSA._run_logged_command(
        `sh -c "exec 3>&1; printf held; sleep 4 >&3 &"`,
        held_stdout_log,
        held_stderr_log)
    @test time() - started < 2.0
    @test read(held_stdout_log, String) == "held"
    @test isempty(read(held_stderr_log, String))

    @test_throws ErrorException Iduna.ThorAxeMSA._run_logged_command(
        `sh -c "exit 2"`,
        joinpath(tmp, "failed_stdout.log"),
        joinpath(tmp, "failed_stderr.log"))

    query_workdir = joinpath(tmp, "query")
    tmp_gene_dir = joinpath(query_workdir, "ENSG00000000001")
    mkpath(joinpath(tmp_gene_dir, "Ensembl"))
    write(joinpath(tmp_gene_dir, "Ensembl", "ensembl_version.csv"), "version\n")
    write(joinpath(tmp_gene_dir, "Ensembl", "tree.nh"), "(a,b);\n")
    missing_outputs = Iduna.ThorAxeMSA._missing_transcript_query_outputs(
        tmp_gene_dir)
    @test "sequences.fasta" in missing_outputs
    @test !("ensembl_version.csv" in missing_outputs)
    @test !("tree.nh" in missing_outputs)

    no_ensembl = joinpath(query_workdir, "NO_ENSEMBL")
    mkpath(no_ensembl)
    @test Iduna.ThorAxeMSA._missing_transcript_query_outputs(no_ensembl) ==
          collect(Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES)
    empty_ensembl = joinpath(query_workdir, "EMPTY_ENSEMBL")
    mkpath(joinpath(empty_ensembl, "Ensembl"))
    @test Iduna.ThorAxeMSA._missing_transcript_query_outputs(empty_ensembl) ==
          collect(Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES)
    complete_ensembl = joinpath(query_workdir, "COMPLETE_ENSEMBL")
    write_test_ensembl_bundle(complete_ensembl)
    @test isempty(Iduna.ThorAxeMSA._missing_transcript_query_outputs(
        complete_ensembl))
end

@testset "ThorAxe environment lifetime" begin
    mktempdir() do tmp
        thoraxe = Iduna.ThorAxeMSA.ThorAxe
        conda = thoraxe.CondaPkg
        # Warm resolution before scheduling concurrent calls. The real wrapper still
        # activates/restores Conda for every call, even with a harmless runner.
        active_env = conda.withenv(() -> copy(ENV))
        original_env = copy(ENV)
        original_frozen = conda.STATE.frozen
        original_dir = pwd()
        target = Iduna.ResolvedTarget(;
            input_id = "ENST", input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG", transcript_id = "ENST")
        input = write_test_ensembl_bundle(joinpath(tmp, "input"))

        function pid_call(workdir, keep, callback)
            actual_thoraxe = (input_dir, output_dir;
                identity, specieslist, phylosofs, runner) -> thoraxe.thoraxe(
                input_dir, output_dir; identity, specieslist, phylosofs,
                aligner = "unused-test-aligner", runner = command -> begin
                    callback(command)
                    write_fake_thoraxe_dir(joinpath(output_dir, "thoraxe"))
                    nothing
                end)
            return Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
                target, input, workdir, 35.0, nothing, keep ? 0 : 1;
                overwrite = true, keep_thoraxe_dir = keep, thoraxe_fn = actual_thoraxe)
        end

        function query_call(workdir, callback; species = nothing)
            mkpath(workdir)
            return Iduna.ThorAxeMSA._run_transcript_query_once(
                "ENSG", workdir, species, nothing,
                joinpath(workdir, "stdout.log"), joinpath(workdir, "stderr.log");
                runner = callback)
        end

        function check_active(_command)
            env_matches = copy(ENV) == active_env
            @test env_matches
            @test conda.STATE.frozen
            return nothing
        end

        pid_call(joinpath(tmp, "warm_kept"), true, check_active)
        pid_call(joinpath(tmp, "warm_temporary"), false, check_active)
        query_call(joinpath(tmp, "warm_query"), check_active)
        query_call(joinpath(tmp, "warm_query_species"), check_active;
            species = "homo_sapiens")
        @test_throws ErrorException query_call(joinpath(tmp, "failed_query"),
            command -> error("expected query runner failure"))
        env_restored = copy(ENV) == original_env
        @test env_restored
        @test conda.STATE.frozen == original_frozen
        @test pwd() == original_dir
        for fail_first in (false, true)
            first_entered = Channel{NamedTuple}(1)
            second_entered = Channel{NamedTuple}(1)
            second_finished = Channel{NamedTuple}(1)
            second_attempting = Channel{Nothing}(1)
            release_first = Channel{Nothing}(1)
            release_second = Channel{Nothing}(1)
            snapshot() = (; env_matches = copy(ENV) == active_env,
                frozen = conda.STATE.frozen, directory = pwd())
            first_callback = command -> begin
                put!(first_entered, snapshot())
                take!(release_first)
                fail_first && error("expected runner failure")
                nothing
            end
            second_callback = command -> begin
                put!(second_entered, snapshot())
                take!(release_second)
                put!(second_finished, snapshot())
                nothing
            end
            first_workdir = joinpath(tmp, "first_$(fail_first)")
            second_workdir = joinpath(tmp, "second_$(fail_first)")
            first_call = Threads.@spawn try
                pid_call(first_workdir, !fail_first, first_callback)
                nothing
            catch err
                err
            end
            second_call = nothing
            try
                @test timedwait(() -> isready(first_entered), 30) == :ok
                isready(first_entered) || error("First ThorAxe callback did not start.")
                first_state = take!(first_entered)
                @test first_state.env_matches
                @test first_state.frozen
                second_call = Threads.@spawn begin
                    put!(second_attempting, nothing)
                    if fail_first
                        pid_call(second_workdir, true, second_callback)
                    else
                        query_call(second_workdir, second_callback)
                    end
                end
                take!(second_attempting)
                # A owns activation until its runner and Conda cleanup both finish.
                @test timedwait(() -> isready(second_entered), 0.5) == :timed_out
                put!(release_first, nothing)
                first_result = fetch(first_call)
                if fail_first
                    @test first_result isa ErrorException
                    @test sprint(showerror, first_result) == "expected runner failure"
                else
                    @test first_result === nothing
                end
                @test timedwait(() -> isready(second_entered), 30) == :ok
                isready(second_entered) || error("Second ThorAxe callback did not start.")
                second_state = take!(second_entered)
                @test second_state.env_matches
                @test second_state.frozen
                env_matches = copy(ENV) == active_env
                @test env_matches
                @test conda.STATE.frozen
                @test pwd() == (fail_first ? original_dir : second_workdir)
            finally
                isready(release_first) || put!(release_first, nothing)
                put!(release_second, nothing)
                wait(first_call)
                second_call === nothing || wait(second_call)
            end
            finished_state = take!(second_finished)
            @test finished_state.env_matches
            @test finished_state.frozen
            env_restored = copy(ENV) == original_env
            @test env_restored
            @test conda.STATE.frozen == original_frozen
            @test pwd() == original_dir
            # A later call must also succeed after an earlier runner throws.
            query_call(joinpath(tmp, "after_$(fail_first)"), check_active)
            env_restored = copy(ENV) == original_env
            @test env_restored
            @test conda.STATE.frozen == original_frozen
            @test pwd() == original_dir
        end
    end
end

mktempdir() do tmp
    source = joinpath(tmp, "cached_input")
    ensembl = joinpath(source, "Ensembl")
    mkpath(ensembl)
    for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
        write(joinpath(ensembl, file), "x\n")
    end
    target = Iduna.ResolvedTarget(;
        input_id = "ENST00000000001.1",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1")
    workdir = joinpath(tmp, "work")
    copied = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
        cached_input_dir = source,
        overwrite = true)
    @test copied == joinpath(workdir, "thoraxe_input")
    @test isfile(joinpath(copied, "Ensembl", "sequences.fasta"))
    @test isdir(source)
    reused = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
        cached_input_dir = source,
        overwrite = false)
    @test reused == copied
    write(joinpath(source, "Ensembl", "sequences.fasta"), "changed\n")
    recopied = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
        cached_input_dir = source,
        overwrite = false)
    @test recopied == copied
    @test read(joinpath(recopied, "Ensembl", "sequences.fasta"), String) == "changed\n"

    @test_throws ErrorException Iduna.ThorAxeMSA._ensure_transcript_query(
        target, joinpath(tmp, "bad_work");
        cached_input_dir = joinpath(tmp, "missing"),
        overwrite = true)

    invalid_metadata_source = joinpath(tmp, "invalid_metadata_source")
    write_test_ensembl_bundle(invalid_metadata_source)
    mkpath(joinpath(invalid_metadata_source,
        Iduna.ThorAxeMSA._TRANSCRIPT_QUERY_METADATA_FILE))
    failed_copy_workdir = joinpath(tmp, "failed_copy_work")
    @test_throws SystemError Iduna.ThorAxeMSA._ensure_transcript_query(
        target, failed_copy_workdir;
        cached_input_dir = invalid_metadata_source,
        overwrite = true)
    failed_input_state = Iduna.Utils._read_stage_state(
        Iduna.ThorAxeMSA._thoraxe_input_stage_dir(failed_copy_workdir))
    @test failed_input_state["status"] == "failed"
    @test failed_input_state["exception"]["type"] == "SystemError"

    direct_workdir = joinpath(tmp, "direct_work")
    direct_input = Iduna.ThorAxeMSA._thoraxe_input_dir(direct_workdir)
    write_test_ensembl_bundle(direct_input)
    direct_metadata = Iduna.ThorAxeMSA._expected_transcript_query_metadata(target;
        specieslist = "homo_sapiens",
        orthology = "1:1",
        source_kind = "transcript_query")
    Iduna.ThorAxeMSA._write_transcript_query_metadata!(
        direct_input, direct_metadata)
    @test Iduna.ThorAxeMSA._ensure_transcript_query(target, direct_workdir;
        specieslist = "homo_sapiens",
        orthology = "1:1") == direct_input

    queried_workdir = joinpath(tmp, "queried_work")
    mkpath(queried_workdir)
    query_calls = Ref(0)
    query_runner = command -> begin
        query_calls[] += 1
        write_test_ensembl_bundle(joinpath(pwd(), "ENSG00000000001"))
        nothing
    end
    query_logs,
    rebuilt_input = Test.collect_test_logs() do
        Iduna.ThorAxeMSA._ensure_transcript_query(
            target, queried_workdir;
            specieslist = "homo_sapiens",
            orthology = "1:1",
            transcript_query_runner = query_runner,
            sleep_fn = seconds -> nothing)
    end
    running_query_log = only([log
                              for log in query_logs
                              if log.message == "Running ThorAxe transcript_query."])
    running_query_kwargs = Dict(running_query_log.kwargs)
    @test haskey(running_query_kwargs, :logs_dir)
    @test !haskey(running_query_kwargs, :stdout_log)
    @test !haskey(running_query_kwargs, :stderr_log)
    @test !haskey(running_query_kwargs, :orthology)
    ready_query_log = only([log
                            for log in query_logs
                            if log.message ==
                               "ThorAxe transcript_query input is ready."])
    @test !haskey(Dict(ready_query_log.kwargs), :input_dir)
    @test rebuilt_input == Iduna.ThorAxeMSA._thoraxe_input_dir(queried_workdir)
    @test query_calls[] == 1

    default_runner_workdir = joinpath(tmp, "default_runner_work")
    mkpath(default_runner_workdir)
    factory_calls = Ref(0)
    runner_calls = Ref(0)
    captured_logs = Tuple{String, String}[]
    default_runner_factory = (stdout_log,
        stderr_log) -> begin
        factory_calls[] += 1
        push!(captured_logs, (stdout_log, stderr_log))
        return command -> begin
            runner_calls[] += 1
            write_test_ensembl_bundle(joinpath(pwd(), "ENSG00000000001"))
            nothing
        end
    end
    default_runner_input = Iduna.ThorAxeMSA._ensure_transcript_query(
        target, default_runner_workdir;
        specieslist = "homo_sapiens",
        orthology = "1:1",
        thoraxe_runner_factory = default_runner_factory,
        sleep_fn = seconds -> nothing)
    @test default_runner_input ==
          Iduna.ThorAxeMSA._thoraxe_input_dir(default_runner_workdir)
    @test factory_calls[] == 1
    @test runner_calls[] == 1
    @test only(captured_logs) == (
        joinpath(default_runner_workdir, "logs", "thoraxe",
            "transcript_query_stdout.log"),
        joinpath(default_runner_workdir, "logs", "thoraxe",
            "transcript_query_stderr.log"))
end
