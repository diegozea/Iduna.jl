@testset "Ensembl homology download retries" begin
    response(status::Integer,
        body::AbstractString = "") = Iduna.Utils.HTTP.Response(status; body)

    attempts = Ref(0)
    data = Iduna.ThorAxeMSA._fetch_ensembl_homology_data(
        "homo_sapiens", "ENSG00000000001.1";
        retries = 3,
        sleep_seconds = 0,
        http_get = (url;
            headers,
            retry,
            status_exception) -> begin
            attempts[] += 1
            @test occursin("/homology/id/homo_sapiens/ENSG00000000001", url)
            @test headers == Iduna.ThorAxeMSA._ENSEMBL_JSON_HEADERS
            @test retry == false
            @test status_exception == false
            attempts[] == 1 ? response(500) : response(200, """{"data": []}""")
        end)
    @test get(data, "data", nothing) !== nothing
    @test attempts[] == 2

    attempts[] = 0
    @test_throws ErrorException Iduna.ThorAxeMSA._fetch_ensembl_homology_data(
        "homo_sapiens", "ENSG00000000001.1";
        retries = 3,
        sleep_seconds = 0,
        http_get = (url; kwargs...) -> begin
            attempts[] += 1
            response(400)
        end)
    @test attempts[] == 1
end

@testset "transcript_query receives orthology and effective species list" begin
    mktempdir() do tmp
        captured = Ref{Cmd}()
        runner = command -> (captured[] = command)
        Iduna.ThorAxeMSA._run_transcript_query_once(
            "ENSG00000000001", tmp, "homo_sapiens",
            "homo_sapiens,mus_musculus",
            joinpath(tmp, "stdout.log"), joinpath(tmp, "stderr.log");
            orthology = "1:n",
            runner = runner)
        parts = captured[].exec
        @test "--orthology" in parts
        @test parts[findfirst(==("--orthology"), parts) + 1] == "1:n"
        @test "--specieslist" in parts
        @test parts[findfirst(==("--specieslist"), parts) + 1] ==
              "homo_sapiens,mus_musculus"

        captured_without_species = Ref{Cmd}()
        Iduna.ThorAxeMSA._run_transcript_query_once(
            "ENSG00000000001", tmp, nothing, nothing,
            joinpath(tmp, "stdout_no_species.log"),
            joinpath(tmp, "stderr_no_species.log");
            orthology = "1:1",
            runner = command -> (captured_without_species[] = command))
        no_species_parts = captured_without_species[].exec
        @test "ENSG00000000001" in no_species_parts
        @test "--orthology" in no_species_parts
    end
end

@testset "transcript_query retry helper branches" begin
    mktempdir() do tmp
        gene_core = "ENSG00000000001"
        stdout_log = joinpath(tmp, "stdout.log")
        stderr_log = joinpath(tmp, "stderr.log")
        tmp_gene_dir = joinpath(tmp, gene_core)

        failed_action = Iduna.ThorAxeMSA._transcript_query_attempt_action!(
            gene_core, tmp, "homo_sapiens", "homo_sapiens,mus_musculus",
            stdout_log, stderr_log, tmp_gene_dir;
            attempt = 1,
            attempts = 1,
            orthology = "1:1",
            runner = command -> nothing)
        @test failed_action === :failed

        invalid_action = @test_logs (:warn, r"invalid bundle") match_mode=:any Iduna.ThorAxeMSA._transcript_query_attempt_action!(
            gene_core, tmp, "homo_sapiens", "homo_sapiens,mus_musculus",
            stdout_log, stderr_log, tmp_gene_dir;
            attempt = 1,
            attempts = 2,
            orthology = "1:1",
            runner = command -> nothing)
        @test invalid_action === :retry

        retry_calls = Ref(0)
        retried_specieslists = String[]
        slept = Float64[]
        retry_runner = command -> begin
            retry_calls[] += 1
            parts = command.exec
            specieslist_index = findfirst(==("--specieslist"), parts)
            push!(retried_specieslists,
                specieslist_index === nothing ? "" : parts[specieslist_index + 1])
            retry_calls[] == 2 &&
                write_test_ensembl_bundle(joinpath(pwd(), gene_core))
            nothing
        end
        @test_logs (:warn, r"invalid bundle") match_mode=:any Iduna.ThorAxeMSA._run_transcript_query_with_retries!(
            tmp_gene_dir, gene_core, tmp, "homo_sapiens",
            "homo_sapiens,mus_musculus", 2, stdout_log, stderr_log;
            orthology = "1:1",
            runner = retry_runner,
            sleep_fn = seconds -> push!(slept, seconds))
        @test retry_calls[] == 2
        @test retried_specieslists ==
              ["homo_sapiens,mus_musculus", "homo_sapiens,mus_musculus"]
        @test slept == [1.0]

        rm(tmp_gene_dir; recursive = true, force = true)
        bundle_script = """
        sleep(0.02)
        gene = $(repr(gene_core))
        files = $(repr(collect(Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES)))
        mkpath(joinpath(gene, "Ensembl"))
        for file in files
            write(joinpath(gene, "Ensembl", file), "x\\n")
        end
        println("done")
        """
        live_retry_runner = _command -> Iduna.ThorAxeMSA._run_logged_command(
            `$(Base.julia_cmd()) --startup-file=no -e $bundle_script`,
            stdout_log,
            stderr_log;
            live_stdout = true)
        tee_stdout_path = joinpath(tmp, "tee_process_stdout.log")
        tee_stderr_path = joinpath(tmp, "tee_process_stderr.log")
        open(tee_stdout_path, "w") do tee_stdout_io
            open(tee_stderr_path, "w") do tee_stderr_io
                redirect_stdout(tee_stdout_io) do
                    redirect_stderr(tee_stderr_io) do
                        Iduna.ThorAxeMSA._run_transcript_query_with_retries!(
                            tmp_gene_dir, gene_core, tmp, "homo_sapiens",
                            "homo_sapiens,mus_musculus", 1, stdout_log, stderr_log;
                            orthology = "1:1",
                            runner = live_retry_runner,
                            sleep_fn = seconds -> nothing)
                    end
                end
            end
        end
        @test Iduna.ThorAxeMSA._has_valid_ensembl_bundle(tmp_gene_dir)
        @test read(stdout_log, String) == "done\n"
        @test isempty(read(tee_stdout_path, String))
        @test occursin("done\n", read(tee_stderr_path, String))

        rm(tmp_gene_dir; recursive = true, force = true)
        failed_calls = Ref(0)
        Iduna.ThorAxeMSA._run_transcript_query_with_retries!(
            tmp_gene_dir, gene_core, tmp, "homo_sapiens",
            "homo_sapiens,mus_musculus", 1, stdout_log, stderr_log;
            orthology = "1:1",
            runner = command -> (failed_calls[] += 1),
            sleep_fn = seconds -> error("failed attempts should not sleep"))
        @test failed_calls[] == 1

        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "$(gene_core).1",
            transcript_id = "ENST00000000001.1",
            species = "homo_sapiens")
        failed_workdir = joinpath(tmp, "failed_query")
        mkpath(failed_workdir)
        try
            Iduna.ThorAxeMSA._ensure_transcript_query(
                target, failed_workdir;
                specieslist = "homo_sapiens,mus_musculus",
                max_retries = 1,
                orthology = "1:1",
                transcript_query_runner = command -> nothing,
                sleep_fn = seconds -> nothing)
            error("expected transcript_query failure")
        catch err
            @test err isa ErrorException
            @test occursin("transcript_query did not create a valid Ensembl bundle",
                sprint(showerror, err))
            @test occursin("smaller curated specieslist", sprint(showerror, err))
        end
    end
end
