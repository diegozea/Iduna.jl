using MIToS.MSA: AnnotatedSequence

@testset "Non-finite sample and reference scores" begin
    mktempdir() do tmp
        input = joinpath(tmp, "input.fasta")
        write(input, ">ref\nACDE\n>other\nACDF\n")
        for value in (NaN, Inf, -Inf)
            score_fn = (reference, sample; kwargs...) -> (; raw_score = value)
            constant_normalization = (score; kwargs...) -> 90.0
            workdir = joinpath(tmp, "raw_$(value)")
            @test_throws ArgumentError Iduna.EPLI.epli_score(input, workdir,
                copy_aligner; score_fn, normalization_fn = constant_normalization,
                sample_count = 1, sample_seed = 1, progress_enabled = false)
            @test !isfile(joinpath(workdir, "summary.csv"))

            reference_score = (reference, sample; label, kwargs...) -> (;
                raw_score = label == "reference_self" ? value : 10.0)
            @test_throws ArgumentError Iduna.EPLI.epli_score(input,
                joinpath(tmp, "reference_$(value)"), copy_aligner;
                score_fn = reference_score,
                normalization_fn = constant_normalization,
                compute_reference_self_score = true, sample_count = 1,
                sample_seed = 1, progress_enabled = false)

            finite_score = (reference, sample; kwargs...) -> (; raw_score = 10.0)
            normalization = (score; kwargs...) -> value
            @test_throws ArgumentError Iduna.EPLI.epli_score(input,
                joinpath(tmp, "normalized_$(value)"), copy_aligner;
                score_fn = finite_score, normalization_fn = normalization,
                sample_count = 1, sample_seed = 1, progress_enabled = false)
        end

        mixed_score = (reference, sample; label, kwargs...) -> (;
            raw_score = endswith(label, "001") ? NaN : 90.0)
        failed_workdir = joinpath(tmp, "mixed_samples")
        error_type = Threads.threadpoolsize() > 1 ? CompositeException : ArgumentError
        @test_throws error_type Iduna.EPLI.epli_score(input, failed_workdir, copy_aligner;
            score_fn = mixed_score, normalization_fn = Iduna.EPLI.no_normalization,
            sample_count = 2, sample_seed = 1, progress_enabled = false)
        @test !isfile(joinpath(failed_workdir, "scores.csv"))
        @test !isfile(joinpath(failed_workdir, "summary.csv"))

        score_fn = (reference, sample; kwargs...) -> (; raw_score = floatmax(Float64))
        workdir = joinpath(tmp, "extreme_finite")
        result = Iduna.EPLI.epli_score(input, workdir, copy_aligner;
            score_fn, normalization_fn = Iduna.EPLI.no_normalization,
            sample_count = 2, sample_seed = 1, progress_enabled = false)
        @test result.epli == floatmax(Float64)
        @test isfile(joinpath(workdir, "scores.csv"))
        @test isfile(joinpath(workdir, "summary.csv"))
    end
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
    @test result.input_source == "fasta"
    @test result.input_fasta == input
    @test result.aligner_args == string(`--example-arg`)
    @test sort(calls) == ["sequence_subset_001", "sequence_subset_002",
        "sequence_subset_003"]
    @test isfile(result.scores_path)
    @test isfile(result.summary_path)

    scores = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(result.scores_path))
    @test scores.epli_component == fill(31.0, 3)
    @test scores.normalization_score == fill(-1.0, 3)
    summary = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(result.summary_path))
    @test only(summary.input_source) == "fasta"
    @test only(summary.input_fasta) == input
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
    sequences = [
        AnnotatedSequence("ref", "AA--AA"),
        AnnotatedSequence("seq2", "A.AA"),
        AnnotatedSequence("seq3", "AAAT")
    ]
    workdir = joinpath(tmp, "mitos_sequences")
    score_fn = (reference_msa, sample_msa; logs_dir = nothing,
        label = "sample") -> begin
        n_sample = Iduna.EPLI.nsequences(
            Iduna.EPLI.read_file(sample_msa, Iduna.EPLI.FASTA))
        return (; raw_score = Float64(n_sample))
    end
    normalization_fn = (score; reference_score = nothing) -> score.raw_score

    result = Iduna.EPLI.epli_score(sequences, workdir, padded_aligner;
        score_fn,
        normalization_fn,
        sample_count = 2,
        sample_fraction = 0.5,
        sample_seed = 23,
        reference_sequence = "seq2",
        progress_enabled = false)

    @test result.epli == 2.0
    @test result.input_source == "mitos_sequences"
    @test ismissing(result.input_fasta)
    full_input = read(joinpath(workdir, "input", "full_sequences.fasta"),
        String)
    @test occursin(">ref\nAAAA\n", full_input)
    @test occursin(">seq2\nAAA\n", full_input)
    @test !occursin('-', full_input)
    @test !occursin('.', full_input)
    sample_input = read(
        joinpath(workdir, "samples", "sequence_subset_001_sequences.fasta"),
        String)
    @test startswith(sample_input, ">seq2\n")
    @test !occursin('-', sample_input)
    @test !occursin('.', sample_input)

    summary = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(result.summary_path))
    @test only(summary.input_source) == "mitos_sequences"
    @test ismissing(only(summary.input_fasta))

    single_result = Iduna.EPLI.epli_score(sequences[1], joinpath(tmp, "single"),
        padded_aligner;
        score_fn,
        normalization_fn,
        sample_count = 1,
        sample_fraction = 1.0,
        sample_seed = 23,
        progress_enabled = false)
    @test single_result.input_source == "mitos_sequence"
    @test ismissing(single_result.input_fasta)
    @test occursin(">ref\nAAAA\n",
        read(joinpath(tmp, "single", "input", "full_sequences.fasta"), String))

    @test_throws ErrorException Iduna.EPLI._sequence_records(AnnotatedSequence[])
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
    aligned_input = joinpath(tmp, "aligned.fasta")
    unaligned_input = joinpath(tmp, "unaligned.fasta")
    write(aligned_input,
        ">ref\nAA--AA\n" *
        ">seq2\nA---AA\n" *
        ">seq3\nA.A-AA\n")
    write(unaligned_input,
        ">ref\nAAAA\n" *
        ">seq2\nAAA\n" *
        ">seq3\nAAAA\n")
    msa = Iduna.EPLI.read_file(aligned_input, Iduna.EPLI.FASTA)
    score_fn = (reference_msa, sample_msa; logs_dir = nothing,
        label = "sample") -> begin
        n_sample = Iduna.EPLI.nsequences(
            Iduna.EPLI.read_file(sample_msa, Iduna.EPLI.FASTA))
        return (; raw_score = Float64(n_sample))
    end
    normalization_fn = (score; reference_score = nothing) -> score.raw_score

    msa_result = Iduna.EPLI.epli_score(msa, joinpath(tmp, "msa"),
        padded_aligner;
        score_fn,
        normalization_fn,
        sample_count = 2,
        sample_fraction = 0.5,
        sample_seed = 19,
        reference_sequence = "seq2",
        progress_enabled = false)
    fasta_result = Iduna.EPLI.epli_score(unaligned_input, joinpath(tmp, "fasta"),
        padded_aligner;
        score_fn,
        normalization_fn,
        sample_count = 2,
        sample_fraction = 0.5,
        sample_seed = 19,
        reference_sequence = "seq2",
        progress_enabled = false)

    @test msa_result.epli == fasta_result.epli
    @test msa_result.input_source == "mitos_msa"
    @test ismissing(msa_result.input_fasta)
    full_input = read(joinpath(tmp, "msa", "input", "full_sequences.fasta"),
        String)
    @test full_input == read(
        joinpath(tmp, "fasta", "input", "full_sequences.fasta"), String)
    @test occursin(">ref\nAAAA\n", full_input)
    @test occursin(">seq2\nAAA\n", full_input)
    @test !occursin('-', full_input)
    @test !occursin('.', full_input)
    sample_input = read(
        joinpath(tmp, "msa", "samples", "sequence_subset_001_sequences.fasta"),
        String)
    @test startswith(sample_input, ">seq2\n")
    @test !occursin('-', sample_input)
    @test !occursin('.', sample_input)

    summary = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(msa_result.summary_path))
    @test only(summary.input_source) == "mitos_msa"
    @test ismissing(only(summary.input_fasta))
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

@testset "Alignment output integrity" begin
    mktempdir() do tmp
        input = joinpath(tmp, "input.fasta")
        output = joinpath(tmp, "alignment.fasta")
        write(input, ">ref full description\nACD\n>other\nAE\n")
        expected = ">other\nA-E\n>ref full description\nACD\n"
        calls = Ref(0)
        fail_after_write = Ref(false)
        aligner = function (source, destination; kwargs...)
            calls[] += 1
            write(destination, fail_after_write[] ? ">other\nA-" : expected)
            fail_after_write[] && error("interrupted alignment")
            return destination
        end
        run_cached = (; overwrite = false) -> Iduna.EPLI._run_aligner(
            aligner, input, output; logs_dir = nothing, run_label = "integrity",
            overwrite, aligner_args = Cmd(String[]))
        @test run_cached() == output
        @test run_cached() == output
        @test calls[] == 1

        # Same-length sequence changes are invisible to existence/size checks.
        for corrupted in (replace(expected, "A-E" => "A-F"), "")
            write(output, corrupted)
            before = calls[]
            @test run_cached() == output
            @test calls[] == before + 1
            @test read(output, String) == expected
        end
        rm(output)
        @test run_cached() == output

        cache_path = Iduna.EPLI._alignment_cache_path(output)
        for change in (
            state -> delete!(state, "schema_version"),
            state -> (state["schema_version"] = 3.0),
            state -> delete!(state, "output_sha256"),
            state -> (state["output_sha256"] = [1, 2]),
            state -> (state["alignment_path"] = "missing.fasta"),
            state -> (state["alignment_path"] = nothing),
            state -> (state["alignment_path"] = "bad\0path"),
            state -> (state["output_fasta"] = "different.fasta"))
            state = Iduna.EPLI._read_alignment_cache(output)
            change(state)
            Iduna.Utils.write_json(cache_path, state)
            before = calls[]
            @test run_cached() == output
            @test calls[] == before + 1
        end
        for malformed in ("{", "[]", "null", "1")
            write(cache_path, malformed)
            before = calls[]
            @test run_cached() == output
            @test calls[] == before + 1
        end

        fail_after_write[] = true
        @test_throws ErrorException run_cached(; overwrite = true)
        @test !isfile(cache_path)
        fail_after_write[] = false
        before = calls[]
        @test run_cached() == output
        @test calls[] == before + 1
        @test read(output, String) == expected
        @test run_cached() == output
        @test calls[] == before + 1
        @test !isfile(string(cache_path, ".tmp"))
    end

    mktempdir() do tmp
        input = joinpath(tmp, "input.fasta")
        output = joinpath(tmp, "alignment.fasta")
        write(input, ">a\nACD\n>b\nAE\n")
        for bad_output in (
            ">a\nACD\n>b\nAE\n", # unequal original row lengths
            ">a\nACD\n>b\nA-F\n", # substituted residue
            ">a\nACD\n>wrong\nA-E\n", # changed ID
            ">a\nACD\n", # missing record
            ">a\nACD\n>a\nA-E\n", # duplicate ID
            ">a\n\n>b\n\n", # zero-column alignment
            ">\nACD\n", # empty ID
            "ACD\n>a\nACD\n", # sequence before first header
            "", ">a")
            aligner = (source, destination; kwargs...) -> write(destination, bad_output)
            returning_aligner = function (source, destination; kwargs...)
                aligner(source, destination; kwargs...)
                return nothing
            end
            @test_throws ErrorException Iduna.EPLI._run_aligner(
                returning_aligner, input, output; logs_dir = nothing,
                run_label = "invalid", overwrite = true, aligner_args = Cmd(String[]))
            @test !isfile(Iduna.EPLI._alignment_cache_path(output))
        end
        missing_aligner = (source, destination; kwargs...) -> nothing
        rm(output; force = true)
        @test_throws ErrorException Iduna.EPLI._run_aligner(
            missing_aligner, input, output; logs_dir = nothing,
            run_label = "missing", overwrite = false, aligner_args = Cmd(String[]))

        alternate = joinpath(tmp, "alternate.fasta")
        calls = Ref(0)
        alternate_aligner = function (source, destination; kwargs...)
            calls[] += 1
            write(alternate, ">b\na.e\n>a\na\ncd\n")
            return (; fasta_path = alternate)
        end
        for _ in 1:2
            @test Iduna.EPLI._run_aligner(alternate_aligner, input, output;
                logs_dir = nothing, run_label = "alternate", overwrite = false,
                aligner_args = Cmd(String[])) == alternate
        end
        @test calls[] == 1
        @test !isfile(output)

        # Compare raw letters before an MSA reader can map unknown residues to X.
        write(input, ">a\nUOX\n>b\nABZ\n")
        replacing_aligner = function (source, destination; kwargs...)
            write(destination, ">a\nXXX\n>b\nABZ\n")
            return destination
        end
        @test_throws ErrorException Iduna.EPLI._run_aligner(
            replacing_aligner, input, output; logs_dir = nothing,
            run_label = "raw_residues", overwrite = true, aligner_args = Cmd(String[]))
        @test !isfile(Iduna.EPLI._alignment_cache_path(output))
    end

    mktempdir() do tmp
        workdir = mkpath(joinpath(tmp, "before"))
        input = joinpath(workdir, "input.fasta")
        output = joinpath(workdir, "alignment.fasta")
        write(input, ">a\nACD\n>b\nACE\n")
        calls = Ref(0)
        aligner = function (source, destination; kwargs...)
            calls[] += 1
            cp(source, destination; force = true)
            return destination
        end
        Iduna.EPLI._run_aligner(aligner, input, output; logs_dir = nothing,
            run_label = "movable", overwrite = false, aligner_args = Cmd(String[]))
        moved = joinpath(tmp, "after")
        mv(workdir, moved)
        moved_output = joinpath(moved, "alignment.fasta")
        @test Iduna.EPLI._run_aligner(
            aligner, joinpath(moved, "input.fasta"), moved_output;
            logs_dir = nothing, run_label = "movable", overwrite = false,
            aligner_args = Cmd(String[])) == moved_output
        @test calls[] == 1
    end
end

@testset "Unreadable cached alignment" begin
    @test_throws ArgumentError Iduna.EPLI._alignment_output_matches("bad\0path", "")
    if !Sys.isunix() || ccall(:geteuid, Cuint, ()) == 0
        @test_skip false # chmod cannot deny reads reliably on these platforms/users.
    else
        mktempdir() do tmp
            input = joinpath(tmp, "input.fasta")
            output = joinpath(tmp, "alignment.fasta")
            cache_path = Iduna.EPLI._alignment_cache_path(output)
            write(input, ">a\nACD\n>b\nACE\n")
            calls = Ref(0)
            fail_rebuild = Ref(false)
            marker_visible = Ref{Union{Nothing, Bool}}(nothing)
            aligner = function (source, destination; kwargs...)
                calls[] += 1
                if fail_rebuild[]
                    marker_visible[] = isfile(cache_path)
                    error("rebuild callback attempted")
                end
                cp(source, destination; force = true)
                return destination
            end
            run_cached = () -> Iduna.EPLI._run_aligner(aligner, input, output;
                logs_dir = nothing, run_label = "unreadable", overwrite = false,
                aligner_args = Cmd(String[]))
            @test run_cached() == output
            @test calls[] == 1
            original_mode = stat(output).mode & 0o777
            chmod(output, 0o000)
            try
                fail_rebuild[] = true
                @test_throws ErrorException run_cached()
                @test calls[] == 2
                @test marker_visible[] === false
                @test !isfile(cache_path)
            finally
                chmod(output, original_mode)
            end
            fail_rebuild[] = false
            @test run_cached() == output
            @test run_cached() == output
            @test calls[] == 3
        end
    end
end
