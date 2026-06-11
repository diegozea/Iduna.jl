@testset "EPLI helpers" begin
    @test !Iduna.Utils._io_is_tty(IOBuffer())
    @test !Iduna.Utils._io_is_tty(IOContext(IOBuffer()))
    withenv("IDUNA_TEST_TRUTHY" => "Yes") do
        @test Iduna.Utils._env_truthy("IDUNA_TEST_TRUTHY")
    end
    withenv("CI" => "false", "GITHUB_ACTIONS" => "false") do
        @test !Iduna.Utils._terminal_progress_enabled(IOBuffer())
    end
    withenv("CI" => "true", "GITHUB_ACTIONS" => "false") do
        @test !Iduna.Utils._terminal_progress_enabled(IOBuffer())
    end

    score = (; raw_score = 12.5)
    normalized = Iduna.EPLI.no_normalization(score)
    @test normalized.normalized_score == 12.5
    @test ismissing(normalized.normalization_score)
end

mktempdir() do tmp
    reference = joinpath(tmp, "reference.fasta")
    sample = joinpath(tmp, "sample.fasta")
    logs_dir = joinpath(tmp, "logs")
    write(reference, ">ref\nAA\n")
    write(sample, ">sample\nAA\n")

    identity_score = Iduna.EPLI.hhsuite_identity_score(reference, sample;
        logs_dir,
        label = "identity")
    @test Iduna.EPLI.comparable_positions_normalization(
        identity_score).normalized_score == 100.0
    @test isfile(joinpath(logs_dir, "identity_hhalign.out"))
    @test Iduna.EPLI._hhsuite_query_indices("A-A", 1, 2) == [1, 0, 2]

    profile_score = Iduna.EPLI.hhsuite_profile_score(reference, reference;
        label = "profile")
    @test Iduna.EPLI.self_reference_normalization(
        profile_score; reference_score = profile_score).normalized_score == 100.0

    @test Iduna.ThorAxeMSA._hhsuite_query_indices("A-A", 1, 2) == [1, 0, 2]
    parsed = Iduna.ThorAxeMSA._parse_hhsuite_query_segment("Q 1 A-A 2")
    @test parsed.indices == [1, 0, 2]
    positions = Int[]
    codes = Char[]
    Iduna.ThorAxeMSA._append_hhsuite_code_line!(positions, codes, "  ||", 3:4,
        [1, 2])
    @test positions == [1, 2]
    @test codes == ['|', '|']
end

mktempdir() do tmp
    input = joinpath(tmp, "input.fasta")
    output = joinpath(tmp, "output.fasta")
    logs_dir = joinpath(tmp, "logs")
    write(input, ">ref\nAA\n")
    captured = Ref{Any}()
    runner = command -> begin
        captured[] = command
        write(output, ">ref\nAA\n")
        return nothing
    end

    path = Iduna.EPLI.prographmsa_aligner(input, output;
        aligner = "/tmp/ProGraphMSA",
        runner,
        logs_dir,
        run_label = "probe",
        aligner_args = `--extra-flag`)

    @test path == output
    @test isfile(joinpath(logs_dir, "probe_prographmsa_stdout.log"))
    @test isfile(joinpath(logs_dir, "probe_prographmsa_stderr.log"))
    command_text = string(captured[])
    @test occursin("--input_order", command_text)
    @test occursin("--fasta", command_text)
    @test occursin("--output", command_text)
    @test occursin("--extra-flag", command_text)
    @test occursin(output, command_text)
    @test occursin(input, command_text)

    plain_output = joinpath(tmp, "plain_output.fasta")
    plain_runner = command -> begin
        captured[] = command
        write(plain_output, ">ref\nAA\n")
        return nothing
    end
    @test Iduna.EPLI.prographmsa_aligner(input, plain_output;
        aligner = "/tmp/ProGraphMSA",
        runner = plain_runner) == plain_output
end
