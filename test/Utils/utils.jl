@testset "Utils" begin
    @test Iduna.Utils.id_kind("P20963") === :uniprot
    @test Iduna.Utils.id_kind("A0A0C4DGY3") === :uniprot
    @test Iduna.Utils.id_kind("A0A0C4DGY3-1") === :uniprot
    @test Iduna.Utils.id_kind("ENST00000362089.10") === :ensembl_transcript
    @test_throws ErrorException Iduna.Utils.id_kind("not-an-id")
    @test Iduna.Utils.format_pid(10) == "10.0"
    @test Iduna.Utils.format_pid(12.34) == "12.3"

    plain = Iduna.Utils.HTTP.Response(200, Vector{UInt8}("plain body"))
    @test Iduna.Utils.decode_body(plain) == "plain body"
    gzip_body = UInt8[
    31, 139, 8, 0, 0, 0, 0, 0, 0, 3, 75, 175, 202, 44, 40, 72, 77, 81, 72,
    202, 79, 169, 4, 0, 196, 225, 85, 59, 12, 0, 0, 0
]
    gzipped = Iduna.Utils.HTTP.Response(200, gzip_body)
    @test Iduna.Utils.decode_body(gzipped) == "gzipped body"

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
