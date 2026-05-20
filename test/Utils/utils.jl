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
    end
end
