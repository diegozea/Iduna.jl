@testset "EPLI" begin
    function copy_aligner(input_fasta, output_fasta; logs_dir = nothing,
            run_label = "run", aligner_args = Cmd(String[]))
        mkpath(dirname(output_fasta))
        cp(input_fasta, output_fasta; force = true)
        return output_fasta
    end

    function padded_aligner(input_fasta, output_fasta; logs_dir = nothing,
            run_label = "run", aligner_args = Cmd(String[]))
        records = Tuple{String, String}[]
        current_name = nothing
        current_seq = String[]
        for line in eachline(input_fasta)
            stripped = strip(line)
            isempty(stripped) && continue
            if startswith(stripped, '>')
                if current_name !== nothing
                    push!(records, (current_name, join(current_seq)))
                end
                current_name = stripped[2:end]
                empty!(current_seq)
            else
                push!(current_seq, stripped)
            end
        end
        if current_name !== nothing
            push!(records, (current_name, join(current_seq)))
        end
        width = maximum(length(seq) for (_name, seq) in records)
        mkpath(dirname(output_fasta))
        open(output_fasta, "w") do io
            for (name, seq) in records
                println(io, ">", name)
                println(io, rpad(seq, width, '-'))
            end
        end
        return output_fasta
    end

    function jll_extension_test_dependencies_available()
        return Base.find_package("MAFFT_jll") !== nothing &&
               Base.find_package("ClustalO_jll") !== nothing
    end

    mktempdir() do tmp
        input = joinpath(tmp, "sequences.fasta")
        write(input,
            ">ref\nAAAA\n" *
            ">seq2\nAAAT\n" *
            ">seq3\nAATT\n" *
            ">seq4\nTTTT\n")
        calls = String[]
        calls_lock = ReentrantLock()
        score_fn = (reference_msa,
            sample_msa;
            logs_dir = nothing,
            label = "sample") -> begin
            lock(calls_lock) do
                push!(calls, String(label))
            end
            n_sample = Iduna.EPLI.nsequences(
                Iduna.EPLI.read_file(sample_msa, Iduna.EPLI.FASTA))
            return (;
                raw_score = 10.0 * n_sample,
                matched_positions = 10 * n_sample,
                comparable_positions = 100)
        end
        normalization_fn = (score;
            reference_score = nothing) -> (;
            normalized_score = score.raw_score + 1,
            normalization_score = -1.0)

        result = Iduna.EPLI.epli_score(input, joinpath(tmp, "run"), copy_aligner;
            score_fn,
            normalization_fn,
            sample_count = 3,
            sample_fraction = 0.5,
            sample_seed = 7,
            aligner_args = `--example-arg`,
            progress_enabled = false)

        @test result.epli == 31.0
        @test result.n_samples == 3
        @test result.sample_seed == UInt64(7)
        @test result.aligner_args == string(`--example-arg`)
        @test sort(calls) == ["sequence_subset_001", "sequence_subset_002",
            "sequence_subset_003"]
        @test isfile(result.scores_path)
        @test isfile(result.summary_path)

        scores = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(result.scores_path))
        @test scores.epli_component == fill(31.0, 3)
        @test scores.normalization_score == fill(-1.0, 3)
        summary = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(result.summary_path))
        @test only(summary.aligner_args) == string(`--example-arg`)

        sample_fasta = joinpath(tmp, "run", "samples", "sequence_subset_001_sequences.fasta")
        @test startswith(read(sample_fasta, String), ">ref\n")

        repeat_result = Iduna.EPLI.epli_score(input, joinpath(tmp, "repeat"), copy_aligner;
            score_fn,
            normalization_fn,
            sample_count = 3,
            sample_fraction = 0.5,
            sample_seed = 7,
            progress_enabled = false)
        for idx in 1:3
            label = "sequence_subset_$(lpad(string(idx), 3, '0'))"
            @test read(joinpath(tmp, "run", "samples", "$(label)_sequences.fasta"), String) ==
                  read(joinpath(tmp, "repeat", "samples", "$(label)_sequences.fasta"),
                String)
        end
        @test repeat_result.epli == result.epli
    end

    mktempdir() do tmp
        input = joinpath(tmp, "sequences.fasta")
        write(input, ">ref\nAAAA\n>seq2\nAAAT\n")
        seen_labels = String[]
        seen_args = String[]
        strict_aligner = (input_fasta, output_fasta; logs_dir = nothing,
            run_label::AbstractString, aligner_args::Cmd) -> begin
            push!(seen_labels, String(run_label))
            push!(seen_args, string(aligner_args))
            return copy_aligner(input_fasta, output_fasta;
                logs_dir, run_label, aligner_args)
        end
        score_fn = (reference_msa, sample_msa; logs_dir = nothing,
            label = "sample") -> (; raw_score = 1.0)
        normalization_fn = (score; reference_score = nothing) -> score.raw_score

        result = Iduna.EPLI.epli_score(input, joinpath(tmp, "strict"),
            strict_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = 17,
            progress_enabled = false)

        @test result.aligner_args == string(Cmd(String[]))
        @test seen_labels == ["full", "sequence_subset_001"]
        @test seen_args == fill(string(Cmd(String[])), 2)
    end

    mktempdir() do tmp
        input = joinpath(tmp, "unaligned.fasta")
        workdir = joinpath(tmp, "unaligned")
        write(input,
            ">ref\nAA--AA\n" *
            ">seq2\nA\n" *
            ">seq3\nA.A\n")
        score_fn = (reference_msa, sample_msa; logs_dir = nothing,
            label = "sample") -> (; raw_score = 5.0)
        normalization_fn = (score; reference_score = nothing) -> score.raw_score

        result = Iduna.EPLI.epli_score(input, workdir, padded_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = 13,
            progress_enabled = false)

        @test result.epli == 5.0
        full_input = read(joinpath(workdir, "input", "full_sequences.fasta"), String)
        sample_input = read(
            joinpath(workdir, "samples",
                "sequence_subset_001_sequences.fasta"), String)
        @test occursin(">ref\nAAAA\n", full_input)
        @test !occursin('-', full_input)
        @test !occursin('.', full_input)
        @test !occursin('-', sample_input)
        @test !occursin('.', sample_input)
    end

    mktempdir() do tmp
        input = joinpath(tmp, "sequences.fasta")
        write(input,
            ">ref\nAAAA\n" *
            ">seq2\nAAAT\n" *
            ">seq3\nAATT\n")
        score_fn = (reference_msa, sample_msa; logs_dir = nothing,
            label = "sample") -> (; raw_score = 7.0)
        scalar_normalization = (score; reference_score = nothing) -> score.raw_score
        namedtuple_aligner = (input_fasta, output_fasta; logs_dir = nothing,
            run_label = "run", aligner_args = Cmd(String[])) -> begin
            mkpath(dirname(output_fasta))
            cp(input_fasta, output_fasta; force = true)
            return (; fasta_path = output_fasta)
        end

        result = Iduna.EPLI.epli_score(input, joinpath(tmp, "named"),
            namedtuple_aligner;
            score_fn,
            normalization_fn = scalar_normalization,
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = 11,
            reference_sequence = "seq2",
            progress_enabled = false)
        @test result.epli == 7.0
        sample_fasta = joinpath(tmp, "named", "samples",
            "sequence_subset_001_sequences.fasta")
        @test startswith(read(sample_fasta, String), ">seq2\n")

        bad_aligner = (input_fasta, output_fasta; logs_dir = nothing,
            run_label = "run", aligner_args = Cmd(String[])) -> 42
        @test_throws ErrorException Iduna.EPLI.epli_score(input, joinpath(tmp, "bad"),
            bad_aligner;
            score_fn,
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = 11,
            progress_enabled = false)
        @test_throws ErrorException Iduna.EPLI.epli_score(input,
            joinpath(tmp, "missing_reference"), copy_aligner;
            score_fn,
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = 11,
            reference_sequence = "missing",
            progress_enabled = false)
        @test_throws ErrorException Iduna.EPLI._normalization_result((value = 1.0,))
        @test_throws ErrorException Iduna.EPLI._normalization_result("bad")
    end

    mktempdir() do tmp
        input = joinpath(tmp, "sequences.fasta")
        workdir = joinpath(tmp, "cached")
        write(input, ">ref\nAAAA\n>seq2\nAAAT\n>seq3\nAATT\n>seq4\nTTTT\n")
        calls = String[]
        cached_aligner = (input_fasta, output_fasta; logs_dir = nothing,
            run_label = "run", aligner_args = Cmd(String[])) -> begin
            push!(calls, String(run_label))
            return copy_aligner(input_fasta, output_fasta;
                logs_dir, run_label, aligner_args)
        end
        score_fn = (reference_msa, sample_msa; logs_dir = nothing,
            label = "sample") -> (; raw_score = 1.0)
        normalization_fn = (score; reference_score = nothing) -> score.raw_score
        Iduna.EPLI.epli_score(input, workdir, cached_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = 3,
            progress_enabled = false)
        @test calls == ["full", "sequence_subset_001"]

        cached = Iduna.EPLI.epli_score(input, workdir, cached_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 1.0,
            sample_seed = 3,
            progress_enabled = false)
        @test cached.epli == 1.0
        @test calls == ["full", "sequence_subset_001"]

        Iduna.EPLI.epli_score(input, workdir, cached_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 0.5,
            sample_seed = 3,
            progress_enabled = false)
        @test calls == ["full", "sequence_subset_001", "sequence_subset_001"]

        arg_result = Iduna.EPLI.epli_score(input, workdir, cached_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 0.5,
            sample_seed = 3,
            aligner_args = `--alpha`,
            progress_enabled = false)
        @test arg_result.aligner_args == string(`--alpha`)
        @test calls == ["full", "sequence_subset_001", "sequence_subset_001",
            "full", "sequence_subset_001"]

        Iduna.EPLI.epli_score(input, workdir, cached_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 0.5,
            sample_seed = 3,
            aligner_args = `--alpha`,
            progress_enabled = false)
        @test calls == ["full", "sequence_subset_001", "sequence_subset_001",
            "full", "sequence_subset_001"]

        Iduna.EPLI.epli_score(input, workdir, cached_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 0.5,
            sample_seed = 3,
            aligner_args = `--beta`,
            progress_enabled = false)
        @test calls == ["full", "sequence_subset_001", "sequence_subset_001",
            "full", "sequence_subset_001", "full", "sequence_subset_001"]

        changed_calls = String[]
        changed_aligner = (input_fasta, output_fasta; logs_dir = nothing,
            run_label = "run", aligner_args = Cmd(String[])) -> begin
            push!(changed_calls, String(run_label))
            return copy_aligner(input_fasta, output_fasta;
                logs_dir, run_label, aligner_args)
        end
        Iduna.EPLI.epli_score(input, workdir, changed_aligner;
            score_fn,
            normalization_fn,
            sample_count = 1,
            sample_fraction = 0.5,
            sample_seed = 3,
            progress_enabled = false)
        @test changed_calls == ["full", "sequence_subset_001"]

        bad_cache_output = joinpath(workdir, "samples", "bad_cache_msa.fasta")
        write(Iduna.EPLI._alignment_cache_path(bad_cache_output), "{")
        @test Iduna.EPLI._read_alignment_cache(bad_cache_output) === nothing
    end

    mktempdir() do tmp
        input = joinpath(tmp, "sequences.fasta")
        write(input, ">ref\nAAAA\n>seq2\nAAAT\n")
        labels = String[]
        labels_lock = ReentrantLock()
        score_fn = (reference_msa, sample_msa; logs_dir = nothing,
            label = "sample") -> begin
            lock(labels_lock) do
                push!(labels, String(label))
            end
            raw = label == "reference_self" ? 20.0 : 10.0
            return (; raw_score = raw)
        end

        result = Iduna.EPLI.epli_score(input, joinpath(tmp, "self"), copy_aligner;
            score_fn,
            normalization_fn = Iduna.EPLI.self_reference_normalization,
            sample_count = 2,
            sample_fraction = 1.0,
            sample_seed = 9,
            progress_enabled = false)

        @test labels[1] == "reference_self"
        @test sort(labels[2:end]) == ["sequence_subset_001", "sequence_subset_002"]
        @test result.epli == 50.0
        scores = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(result.scores_path))
        @test scores.normalization_score == fill(20.0, 2)
    end

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

    @testset "JLL aligner extensions" begin
        @test isdefined(Iduna.EPLI, :mafft_aligner)
        @test isdefined(Iduna.EPLI, :clustalo_aligner)
        @test !(:mafft_aligner in names(Iduna.EPLI))
        @test !(:clustalo_aligner in names(Iduna.EPLI))

        if jll_extension_test_dependencies_available()
            @test Base.get_extension(Iduna, :IdunaMAFFTExt) === nothing
            @eval using MAFFT_jll
            @test Base.get_extension(Iduna, :IdunaMAFFTExt) !== nothing

            @test Base.get_extension(Iduna, :IdunaClustalOExt) === nothing
            @eval using ClustalO_jll
            @test Base.get_extension(Iduna, :IdunaClustalOExt) !== nothing

            mktempdir() do tmp
                input = joinpath(tmp, "input.fasta")
                output = joinpath(tmp, "mafft_output.fasta")
                logs_dir = joinpath(tmp, "logs")
                write(input, ">ref\nAA\n")
                captured = Ref{Any}()
                runner = command -> begin
                    captured[] = command
                    write(output, ">ref\nAA\n")
                    return nothing
                end

                path = Iduna.EPLI.mafft_aligner(input, output;
                    aligner = "/tmp/mafft",
                    runner,
                    logs_dir,
                    run_label = "probe",
                    aligner_args = `--maxiterate 1000 --localpair`)

                @test path == output
                @test isfile(output)
                @test isfile(joinpath(logs_dir, "probe_mafft_stderr.log"))
                command_text = string(captured[])
                @test occursin("--quiet", command_text)
                @test occursin("--inputorder", command_text)
                @test occursin("--maxiterate", command_text)
                @test occursin("1000", command_text)
                @test occursin("--localpair", command_text)
                @test occursin(input, command_text)
                @test occursin(output, command_text)

                plain_output = joinpath(tmp, "mafft_plain_output.fasta")
                plain_runner = command -> begin
                    captured[] = command
                    write(plain_output, ">ref\nAA\n")
                    return nothing
                end
                @test Iduna.EPLI.mafft_aligner(input, plain_output;
                    aligner = "/tmp/mafft",
                    runner = plain_runner) == plain_output
            end

            mktempdir() do tmp
                input = joinpath(tmp, "input.fasta")
                output = joinpath(tmp, "clustalo_output.fasta")
                logs_dir = joinpath(tmp, "logs")
                write(input, ">ref\nAA\n")
                captured = Ref{Any}()
                runner = command -> begin
                    captured[] = command
                    write(output, ">ref\nAA\n")
                    return nothing
                end

                path = Iduna.EPLI.clustalo_aligner(input, output;
                    aligner = "/tmp/clustalo",
                    runner,
                    logs_dir,
                    run_label = "probe",
                    aligner_args = `--threads=4 --auto`)

                @test path == output
                @test isfile(output)
                @test isfile(joinpath(logs_dir, "probe_clustalo_stdout.log"))
                @test isfile(joinpath(logs_dir, "probe_clustalo_stderr.log"))
                command_text = string(captured[])
                @test occursin("--infile", command_text)
                @test occursin("--outfile", command_text)
                @test occursin("--outfmt=fasta", command_text)
                @test occursin("--force", command_text)
                @test occursin("--output-order=input-order", command_text)
                @test occursin("--threads=4", command_text)
                @test occursin("--auto", command_text)
                @test occursin(input, command_text)
                @test occursin(output, command_text)

                plain_output = joinpath(tmp, "clustalo_plain_output.fasta")
                plain_runner = command -> begin
                    captured[] = command
                    write(plain_output, ">ref\nAA\n")
                    return nothing
                end
                @test Iduna.EPLI.clustalo_aligner(input, plain_output;
                    aligner = "/tmp/clustalo",
                    runner = plain_runner) == plain_output
            end

            mktempdir() do tmp
                input = joinpath(tmp, "sequences.fasta")
                write(input,
                    ">ref\nAAAA\n" *
                    ">seq2\nAAAT\n" *
                    ">seq3\nAATT\n")
                mafft_copy_aligner = (input_fasta, output_fasta;
                    logs_dir = nothing,
                    run_label = "run",
                    aligner_args = Cmd(String[])) -> begin
                    runner = command -> begin
                        padded_aligner(input_fasta, output_fasta;
                            logs_dir, run_label, aligner_args)
                        return nothing
                    end
                    return Iduna.EPLI.mafft_aligner(input_fasta, output_fasta;
                        logs_dir,
                        run_label,
                        aligner_args,
                        runner,
                        aligner = "/tmp/mafft")
                end
                score_fn = (reference_msa, sample_msa; logs_dir = nothing,
                    label = "sample") -> (; raw_score = 3.0)

                result = Iduna.EPLI.epli_score(input, joinpath(tmp, "mafft_epli"),
                    mafft_copy_aligner;
                    score_fn,
                    normalization_fn = Iduna.EPLI.no_normalization,
                    sample_count = 1,
                    sample_fraction = 1.0,
                    sample_seed = 5,
                    aligner_args = `--thread 1`,
                    progress_enabled = false)

                @test isfinite(result.epli)
                @test isfile(result.scores_path)
                @test isfile(result.summary_path)
                @test result.aligner_args == string(`--thread 1`)
            end
        else
            @info "Skipping optional JLL aligner extension tests; MAFFT_jll and ClustalO_jll are not available in this environment."
        end
    end
end
