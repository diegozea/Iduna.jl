@testset "Utils" begin
    @test Iduna.Utils.id_kind("P20963") === :uniprot
    @test Iduna.Utils.id_kind("A0A0C4DGY3") === :uniprot
    @test Iduna.Utils.id_kind("A0A0C4DGY3-1") === :uniprot
    @test Iduna.Utils.id_kind("ENST00000362089.10") === :ensembl_transcript
    @test_throws ErrorException Iduna.Utils.id_kind("not-an-id")
    @test Iduna.Utils.format_pid(10) == "10.0"
    @test Iduna.Utils.format_pid(12.34) == "12.3"
    @test Iduna.Utils.format_pid_dir(10) == "pid_10.00"
    @test Iduna.Utils.format_pid_dir(12.34) == "pid_12.34"

    @testset "s-exon provenance annotations" begin
        mktempdir() do tmp
            sto = joinpath(tmp, "annotated.sto")
            write(sto, """
            # STOCKHOLM 1.0
            #=GF SExonCodeMap "0"=>"1_0","β"=>"2_0","="=>"3_0"
            seq1 ACDE
            seq2 A-DE
            #=GC SExonCode 00β=
            //
            """)
            msa = Iduna.MSAExpansion.read_file(sto, Iduna.MSAExpansion.Stockholm;
                keepinserts = true)
            @test Iduna.Utils.s_exon_codes(msa) == "00β="
            @test Iduna.Utils.s_exon_code_map(msa) ==
                  Dict('0' => "1_0", 'β' => "2_0", '=' => "3_0")
            @test Iduna.Utils.s_exon_codes(msa[:, [1, 3, 4]]) == "0β="

            blocks = joinpath(tmp, "blocks.tsv")
            Iduna.Utils.write_s_exon_blocks_tsv(blocks, msa;
                alignment = "seed",
                pid = 10.0)
            @test read(blocks, String) ==
                  "alignment\tpid\tcode\ts_exon_id\tstart_col\tend_col\tn_columns\n" *
                  "seed\t10.0\t0\t1_0\t1\t2\t2\n" *
                  "seed\t10.0\tβ\t2_0\t3\t3\t1\n" *
                  "seed\t10.0\t=\t3_0\t4\t4\t1\n"
        end
    end

    plain = Iduna.Utils.HTTP.Response(200, Vector{UInt8}("plain body"))
    @test Iduna.Utils.decode_body(plain) == "plain body"
    gzip_body = UInt8[
    31, 139, 8, 0, 0, 0, 0, 0, 0, 3, 75, 175, 202, 44, 40, 72, 77, 81, 72,
    202, 79, 169, 4, 0, 196, 225, 85, 59, 12, 0, 0, 0
]
    gzipped = Iduna.Utils.HTTP.Response(200, gzip_body)
    @test Iduna.Utils.decode_body(gzipped) == "gzipped body"

    @testset "HTTP download retries" begin
        response(status::Integer,
            body::AbstractString = "") = Iduna.Utils.HTTP.Response(status, Vector{UInt8}(body))

        attempts = Ref(0)
        transient_then_success = (url; headers, retry,
            status_exception) -> begin
            attempts[] += 1
            @test url == "https://example.test/data"
            @test headers == ["Accept" => "text/plain"]
            @test retry == false
            @test status_exception == false
            attempts[] == 1 ? response(503) : response(200, "ok")
        end
        resp = Iduna.Utils._http_get_with_retries(
            "https://example.test/data", ["Accept" => "text/plain"];
            retries = 3,
            sleep_seconds = 0,
            http_get = transient_then_success)
        @test resp.status == 200
        @test Iduna.Utils.decode_body(resp) == "ok"
        @test attempts[] == 2

        attempts[] = 0
        non_transient = (url; kwargs...) -> begin
            attempts[] += 1
            response(404)
        end
        @test Iduna.Utils._http_get_with_retries(
            "https://example.test/missing", [];
            retries = 3,
            sleep_seconds = 0,
            http_get = non_transient).status == 404
        @test attempts[] == 1

        attempts[] = 0
        exhausted = (url; kwargs...) -> begin
            attempts[] += 1
            response(503)
        end
        @test Iduna.Utils._http_get_with_retries(
            "https://example.test/unavailable", [];
            retries = 2,
            sleep_seconds = 0,
            http_get = exhausted).status == 503
        @test attempts[] == 2

        attempts[] = 0
        @test Iduna.Utils._http_get_with_retries(
            "https://example.test/once", [];
            retries = 1,
            sleep_seconds = 0,
            http_get = exhausted).status == 503
        @test attempts[] == 1

        throws_unretryable = (url; kwargs...) -> error("network boom")
        @test_throws ErrorException Iduna.Utils._http_get_with_retries(
            "https://example.test/error", [];
            retries = 2,
            sleep_seconds = 0,
            http_get = throws_unretryable)
    end

    identical = Iduna.Utils.protein_alignment_stats("ACDE", "ACDE")
    @test identical.identical
    @test identical.mismatches == 0
    @test identical.insertions == 0
    @test identical.deletions == 0

    mismatch = Iduna.Utils.protein_alignment_stats("ACDF", "ACDE")
    @test !mismatch.identical
    @test mismatch.mismatches == 1

    indel = Iduna.Utils.protein_alignment_stats("ACDEF", "ACDE")
    @test !indel.identical
    @test indel.insertions + indel.deletions == 1

    mktempdir() do tmp
        sto = joinpath(tmp, "names.sto")
        write(sto, "# STOCKHOLM 1.0\n1|ENSG000001.13 AC\nENST00000362089 AC\n//\n")
        msa = Iduna.ResultsValidation.load_msa(sto)
        @test Iduna.Utils.resolve_sequence_name(msa, "ENSG000001.13") ==
              "1|ENSG000001.13"
        @test Iduna.Utils.resolve_sequence_name(msa, "ENSG000001") ==
              "1|ENSG000001.13"
        @test Iduna.Utils.resolve_sequence_name(msa, "ENST00000362089.10") ==
              "ENST00000362089"
        @test Iduna.Utils.resolve_sequence_name(msa, "missing") === nothing
        @test Iduna.Utils.resolve_sequence_name(msa, "missing"; fallback = true) ==
              "1|ENSG000001.13"

        workdir = Iduna.Utils.prepare_output_dir("P20963"; output_dir = joinpath(tmp, "P20963"))
        @test isdir(workdir)
        @test Iduna.Utils.prepare_output_dir(
            "P20963"; workdir, output_dir = workdir) == workdir
        @test_throws ErrorException Iduna.Utils.prepare_output_dir(
            "P20963"; workdir, output_dir = joinpath(tmp, "other"))
        file_output = joinpath(tmp, "file-output")
        write(file_output, "already here")
        @test_throws ErrorException Iduna.Utils.prepare_output_dir(
            "P20963"; output_dir = file_output)

        owned = joinpath(workdir, "owned")
        mkpath(owned)
        @test isdir(owned)
        Iduna.Utils.safe_rm(owned, workdir)
        @test !isdir(owned)
        @test_throws ErrorException Iduna.Utils.safe_rm(tmp, workdir)

        inside_path = joinpath(workdir, "nested", "artifact.txt")
        outside_path = joinpath(tmp, "outside.txt")
        relative_path = joinpath("nested", "artifact.txt")
        @test Iduna.Utils._relative_artifact_path(inside_path, workdir) == relative_path
        @test Iduna.Utils._relative_artifact_path(outside_path, workdir) == outside_path
        @test Iduna.Utils._relative_artifact_path(relative_path, workdir) == relative_path
        @test Iduna.Utils._relative_artifact_path(nothing, workdir) === nothing
        @test Iduna.Utils._resolve_artifact_path(relative_path, workdir) == inside_path
        @test Iduna.Utils._resolve_artifact_path(outside_path, workdir) == outside_path
        @test Iduna.Utils._resolve_artifact_path(nothing, workdir) === nothing

        formatted = Iduna.Utils.format_fasta("seq", "acdefghijklmnopqrstuvwxyz"^3)
        @test startswith(formatted, ">seq\nACDEFGHIJKLMNOPQRSTUVWXYZ")
        @test occursin("\nLMNOPQRSTUVWXYZ\n", formatted)

        fasta_path = joinpath(tmp, "nested", "records.fasta")
        @test Iduna.Utils.write_fasta(
            fasta_path, [("seq1", "acde"), ("seq2", "fghi")]) == fasta_path
        @test read(fasta_path, String) == ">seq1\nACDE\n\n>seq2\nFGHI\n\n"

        text_path = joinpath(tmp, "text", "note.txt")
        @test Iduna.Utils.write_text(text_path, "hello") == text_path
        @test read(text_path, String) == "hello"

        stdout_path = joinpath(tmp, "logs", "stdout.log")
        stderr_path = joinpath(tmp, "logs", "stderr.log")
        Iduna.Utils.run_logged(
            `sh -c "printf out; printf err >&2"`;
            stdout_path,
            stderr_path
        )
        @test read(stdout_path, String) == "out"
        @test read(stderr_path, String) == "err"

        cd_stdout = joinpath(tmp, "logs", "pwd_stdout.log")
        cd_stderr = joinpath(tmp, "logs", "pwd_stderr.log")
        Iduna.Utils.run_logged(`pwd`; stdout_path = cd_stdout, stderr_path = cd_stderr,
            workdir)
        @test realpath(chomp(read(cd_stdout, String))) == realpath(workdir)
        @test isempty(read(cd_stderr, String))

        stage_dir = Iduna.Utils._pipeline_stage_dir(workdir, "example:stage")
        stage_identity = (; input = "A", option = 1)
        stage_output = joinpath(workdir, "stage", "artifact.txt")
        stage_outputs = (; artifact = stage_output)
        missing_stage = Iduna.Utils._classify_stage_state(
            stage_dir, stage_identity, stage_outputs; stage_label = "example")
        @test missing_stage.status === :missing
        @test missing_stage.reusable === false

        mkpath(dirname(stage_output))
        write(stage_output, "artifact")
        stale_stage = Iduna.Utils._classify_stage_state(
            stage_dir, stage_identity, stage_outputs; stage_label = "example")
        @test stale_stage.status === :stale
        @test occursin("no stage_state.json", stale_stage.warning)

        Iduna.Utils._write_stage_state(stage_dir;
            stage = "example",
            stage_key = "example:stage",
            status = :running,
            identity = stage_identity,
            outputs = stage_outputs,
            action = :run,
            workdir)
        running_stage = Iduna.Utils._classify_stage_state(
            stage_dir, stage_identity, stage_outputs; stage_label = "example")
        @test running_stage.status === :unfinished

        Iduna.Utils._write_stage_state(stage_dir;
            stage = "example",
            stage_key = "example:stage",
            status = :done,
            identity = stage_identity,
            outputs = stage_outputs,
            action = :run,
            workdir)
        done_stage = Iduna.Utils._classify_stage_state(
            stage_dir, stage_identity, stage_outputs; stage_label = "example")
        @test done_stage.reusable === true
        @test done_stage.status === :done
        changed_stage = Iduna.Utils._classify_stage_state(
            stage_dir, (; input = "B", option = 1), stage_outputs;
            stage_label = "example")
        @test changed_stage.status === :stale

        summaries = Iduna.Utils.collect_stage_summaries(workdir)
        summary_idx = findfirst(
            summary -> summary.stage_key == "example:stage", summaries)
        @test summary_idx !== nothing
        @test summaries[summary_idx].action == "run"
        @test summaries[summary_idx].outputs["artifact"] ==
              joinpath("stage", "artifact.txt")
        filtered_summaries = Iduna.Utils.collect_stage_summaries(workdir;
            stage_keys = ["example:stage"])
        @test length(filtered_summaries) == 1
        @test only(filtered_summaries).stage_key == "example:stage"
        @test isempty(Iduna.Utils.collect_stage_summaries(workdir;
            stage_keys = ["missing:stage"]))
        @test Iduna.Utils._stage_state_unreadable_message(nothing) ==
              "state file disappeared while reading"
        invalid_stage_dir = joinpath(workdir, "invalid-stage")
        mkpath(invalid_stage_dir)
        write(joinpath(invalid_stage_dir, "stage_state.json"), "{bad json")
        unreadable_state = Iduna.Utils._read_stage_state(invalid_stage_dir)
        @test unreadable_state isa NamedTuple
        @test Iduna.Utils._stage_state_unreadable_message(unreadable_state) !== nothing
        invalid_state = joinpath(workdir, "invalid_stage_state.json")
        write(invalid_state, "{bad json")
        @test Iduna.Utils._stage_summary_from_state_path(invalid_state, workdir) === nothing
    end
end
