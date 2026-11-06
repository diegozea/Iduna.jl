using MIToS.MSA: AnnotatedSequence

@testset "Typed callback contract" begin
    E = Iduna.EPLI
    mktempdir() do tmp
        input = write_epli_test_input(tmp)
        valid_score = constant_score_callback(10)
        function tuple_reference_score(args...; label, kwargs...)
            label == "reference_self" && return (; raw_score = 10)
            return E.RawScore(10)
        end
        cases = (
            ("score_tuple", (; score_fn = (args...; kwargs...) -> (; raw_score = 10))),
            ("reference_tuple",
                (; score_fn = tuple_reference_score,
                    compute_reference_self_score = true)),
            ("normalization_scalar",
                (; score_fn = valid_score,
                    normalization_fn = (score; kwargs...) -> 10.0)),
            ("normalization_tuple",
                (; score_fn = valid_score,
                    normalization_fn = (score; kwargs...) -> (; normalized_score = 10.0)))
        )
        for (label, callbacks) in cases
            workdir = joinpath(tmp, label)
            @test_throws TypeError E.epli_score(input, workdir, copy_aligner;
                callbacks..., sample_count = 1, sample_seed = 1, progress_enabled = false)
            @test !isfile(joinpath(workdir, "scores.csv"))
            @test !isfile(joinpath(workdir, "summary.csv"))
        end
    end
end

@testset "Custom score implementations" begin
    E = Iduna.EPLI
    @test isabstracttype(E.ScoreResult)
    raw = CustomRawScore(Float32(-2))
    @test E.identity_counts(raw) === nothing
    @test E.score_diagnostics(raw) == (;)
    @test_throws ArgumentError E.comparable_positions_normalization(raw)
    # Custom numeric storage must not introduce integer or Float32 overflow.
    for value in (typemax(Int), floatmax(Float32), 3 // 4)
        normalized = E.self_reference_normalization(CustomRawScore(value);
            reference_score = CustomRawScore(value))
        @test normalized.normalized_score ≈ 100.0
        @test normalized.normalization_score === Float64(value)
    end
    for invalid in (NaN, Inf, -Inf)
        @test_throws ArgumentError E.self_reference_normalization(CustomRawScore(invalid);
            reference_score = CustomRawScore(1))
        @test_throws ArgumentError E.self_reference_normalization(CustomRawScore(1);
            reference_score = CustomRawScore(invalid))
    end
    mktempdir() do tmp
        input = write_epli_test_input(tmp)
        score_fn = (args...; label, kwargs...) -> label == "reference_self" ?
                                                  CustomRawScore(4) : raw
        for (normalizer, expected) in ((E.no_normalization, -2.0),
            (E.self_reference_normalization, -50.0))
            result = E.epli_score(input, joinpath(tmp, string(normalizer)), copy_aligner;
                score_fn, normalization_fn = normalizer, sample_count = 1,
                sample_seed = 1, progress_enabled = false)
            @test Iduna.epli_value(result) === expected
            @test result.samples_without_identity_counts == 1
            @test ismissing(only(result.rows).e_value)
        end

        for counts in (E.IdentityCounts(3, 4), E.IdentityCounts(0, 0))
            score = CustomCountScore(-2.0, counts, (; e_value = Inf, note = :preserved))
            function normalize_custom(value; reference_score)
                @test value === score # The normalizer receives the user's own type.
                @test E.score_diagnostics(value).note === :preserved
                return E.comparable_positions_normalization(value)
            end
            result = E.epli_score(
                input, joinpath(tmp, "counts_$(counts.comparable_positions)"),
                copy_aligner; score_fn = (args...; kwargs...) -> score,
                normalization_fn = normalize_custom, sample_count = 1,
                sample_seed = 1, progress_enabled = false)
            row = only(result.rows)
            @test Iduna.epli_value(result) ==
                  (counts.comparable_positions == 0 ? 0.0 : 75.0)
            @test row.raw_score === -2.0
            @test row.matched_positions == counts.matched_positions
            @test row.comparable_positions == counts.comparable_positions
            @test row.no_comparable_positions == iszero(counts.comparable_positions)
            @test row.e_value === Inf
        end
    end
end

@testset "Scoring failures do not produce results" begin
    E = Iduna.EPLI
    mktempdir() do tmp
        input = write_epli_test_input(tmp)
        normalization_calls = Ref(0)
        observed_normalization = (score; kwargs...) -> begin
            normalization_calls[] += 1
            return E.NormalizationResult(90)
        end
        function invalid_reference_score(args...; label, kwargs...)
            value = label == "reference_self" ? NaN : 10
            return E.RawScore(value)
        end
        invalid_normalization = (score; kwargs...) -> E.NormalizationResult(NaN)
        # The constructor tests cover NaN and both infinities. Here each boundary
        # must stop the pipeline before aggregation or output, using NaN as a probe.
        cases = (
            ("sample",
                (; score_fn = constant_score_callback(NaN),
                    normalization_fn = observed_normalization)),
            ("reference",
                (; score_fn = invalid_reference_score,
                    normalization_fn = observed_normalization,
                    compute_reference_self_score = true)),
            ("normalization",
                (; score_fn = constant_score_callback(10),
                    normalization_fn = invalid_normalization))
        )
        for (label, callbacks) in cases
            workdir = joinpath(tmp, label)
            @test_throws ArgumentError E.epli_score(input, workdir, copy_aligner;
                callbacks..., sample_count = 1, sample_seed = 1, progress_enabled = false)
            @test !isfile(joinpath(workdir, "scores.csv"))
            @test !isfile(joinpath(workdir, "summary.csv"))
        end

        # Subtyping alone cannot guarantee valid values. Both callback boundaries
        # must reject bad extension data even when the normalizer would ignore it.
        invalid_custom_scores = (
            ("nan", CustomRawScore(NaN), ArgumentError),
            ("inf", CustomRawScore(Inf), ArgumentError),
            ("negative_inf", CustomRawScore(-Inf), ArgumentError),
            ("overflow", CustomRawScore(big"1e1000"), ArgumentError),
            ("nonnumeric", CustomRawScore("bad"), TypeError),
            ("counts", CustomCountScore(1.0, (; matched = 3, comparable = 2), (;)),
                TypeError),
            ("missing_counts", CustomCountScore(1.0, missing, (;)), TypeError))
        for (name, bad_score, error_type) in invalid_custom_scores
            for (stage, self_score) in (("sample", false), ("reference", true))
                score_fn = (args...; label, kwargs...) -> begin
                    if self_score && label != "reference_self"
                        return E.RawScore(10)
                    end
                    return bad_score
                end
                workdir = joinpath(tmp, "$(stage)_$(name)")
                @test_throws error_type E.epli_score(input, workdir, copy_aligner;
                    score_fn, normalization_fn = observed_normalization,
                    compute_reference_self_score = self_score, sample_count = 1,
                    sample_seed = 1, progress_enabled = false)
                @test !isfile(joinpath(workdir, "scores.csv"))
                @test !isfile(joinpath(workdir, "summary.csv"))
            end
        end
        @test normalization_calls[] == 0

        mixed_score = (args...; label, kwargs...) -> begin
            value = endswith(label, "001") ? NaN : 90
            return E.RawScore(value)
        end
        failed_workdir = joinpath(tmp, "mixed_samples")
        error_type = Threads.threadpoolsize() > 1 ? CompositeException : ArgumentError
        @test_throws error_type E.epli_score(input, failed_workdir, copy_aligner;
            score_fn = mixed_score, normalization_fn = E.no_normalization,
            sample_count = 2, sample_seed = 1, progress_enabled = false)
        @test !isfile(joinpath(failed_workdir, "scores.csv"))
        @test !isfile(joinpath(failed_workdir, "summary.csv"))

        result = E.epli_score(
            input, joinpath(tmp, "extreme_finite"), copy_aligner;
            score_fn = constant_score_callback(floatmax(Float64)),
            normalization_fn = E.no_normalization,
            sample_count = 2, sample_seed = 1, progress_enabled = false)
        @test result.epli.value == floatmax(Float64)
        @test isfile(result.scores_path)
        @test isfile(result.summary_path)
    end
end

@testset "Optional score metadata at callback boundaries" begin
    mktempdir() do tmp
        input = write_epli_test_input(tmp)
        score = Iduna.EPLI.RawScore(Float32(-2);
            diagnostics = (;
                hhalign_score = NaN, probability = missing, e_value = Inf,
                matched_positions = 12, comparable_positions = 0,
                custom_metadata = :preserved))
        score_fn = (reference, sample; kwargs...) -> score
        normalization_fn = (value; kwargs...) -> begin
            @test Iduna.EPLI.score_diagnostics(value).custom_metadata === :preserved
            @test Iduna.EPLI.score_diagnostics(value).matched_positions == 12
            return Iduna.EPLI.NormalizationResult(150)
        end
        result = Iduna.EPLI.epli_score(input, joinpath(tmp, "run"), copy_aligner;
            score_fn, normalization_fn, sample_count = 1, sample_seed = 1,
            progress_enabled = false)
        row = only(result.rows)
        @test row.raw_score === -2.0
        @test row.normalized_score === 150.0
        @test ismissing(row.normalization_score)
        @test ismissing(row.matched_positions)
        @test ismissing(row.comparable_positions)
        @test isnan(row.hhalign_score)
        @test ismissing(row.probability)
        @test row.e_value === Inf
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
        return Iduna.EPLI.IdentityScore(Iduna.EPLI.IdentityCounts(10 * n_sample, 100))
    end
    normalization_fn = (score; reference_score = nothing) -> Iduna.EPLI.NormalizationResult(
        Iduna.EPLI.raw_score(score) + 1, -1)

    result = Iduna.EPLI.epli_score(input, joinpath(tmp, "run"), copy_aligner;
        score_fn,
        normalization_fn,
        sample_count = 3,
        sample_fraction = 0.5,
        sample_seed = 7,
        aligner_args = `--example-arg`,
        progress_enabled = false)

    @test Iduna.epli_value(result) === 31.0
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
    @test repeat_result.epli.value == result.epli.value
end

mktempdir() do tmp
    sequences = [
        AnnotatedSequence("ref", "AA--AA"),
        AnnotatedSequence("seq2", "A.AA"),
        AnnotatedSequence("seq3", "AAAT")
    ]
    workdir = joinpath(tmp, "mitos_sequences")
    score_fn = sample_size_score
    normalization_fn = Iduna.EPLI.no_normalization

    result = Iduna.EPLI.epli_score(sequences, workdir, padded_aligner;
        score_fn,
        normalization_fn,
        sample_count = 2,
        sample_fraction = 0.5,
        sample_seed = 23,
        reference_sequence = "seq2",
        progress_enabled = false)

    @test result.epli.value == 2.0
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
    score_fn = constant_score_callback(1.0)
    normalization_fn = Iduna.EPLI.no_normalization

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
    score_fn = constant_score_callback(5.0)
    normalization_fn = Iduna.EPLI.no_normalization

    result = Iduna.EPLI.epli_score(input, workdir, padded_aligner;
        score_fn,
        normalization_fn,
        sample_count = 1,
        sample_fraction = 1.0,
        sample_seed = 13,
        progress_enabled = false)

    @test result.epli.value == 5.0
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
    score_fn = sample_size_score
    normalization_fn = Iduna.EPLI.no_normalization

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

    @test msa_result.epli.value == fasta_result.epli.value
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
    score_fn = constant_score_callback(7.0)
    normalization_fn = Iduna.EPLI.no_normalization
    namedtuple_aligner = (input_fasta, output_fasta; logs_dir = nothing,
        run_label = "run", aligner_args = Cmd(String[])) -> begin
        mkpath(dirname(output_fasta))
        cp(input_fasta, output_fasta; force = true)
        return (; fasta_path = output_fasta)
    end

    result = Iduna.EPLI.epli_score(input, joinpath(tmp, "named"),
        namedtuple_aligner;
        score_fn,
        normalization_fn,
        sample_count = 1,
        sample_fraction = 1.0,
        sample_seed = 11,
        reference_sequence = "seq2",
        progress_enabled = false)
    @test result.epli.value == 7.0
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
    score_fn = constant_score_callback(1.0)
    normalization_fn = Iduna.EPLI.no_normalization
    function run_cached(; aligner_fn = cached_aligner, sample_fraction = 1.0, kwargs...)
        Iduna.EPLI.epli_score(input, workdir, aligner_fn;
            score_fn, normalization_fn, sample_count = 1, sample_seed = 3,
            sample_fraction, progress_enabled = false, kwargs...)
    end
    run_cached()
    @test calls == ["full", "sequence_subset_001"]

    cached = run_cached()
    @test cached.epli.value == 1.0
    @test calls == ["full", "sequence_subset_001"]

    run_cached(sample_fraction = 0.5)
    @test calls == ["full", "sequence_subset_001", "sequence_subset_001"]

    arg_result = run_cached(sample_fraction = 0.5, aligner_args = `--alpha`)
    @test arg_result.aligner_args == string(`--alpha`)
    @test calls == ["full", "sequence_subset_001", "sequence_subset_001",
        "full", "sequence_subset_001"]

    run_cached(sample_fraction = 0.5, aligner_args = `--alpha`)
    @test calls == ["full", "sequence_subset_001", "sequence_subset_001",
        "full", "sequence_subset_001"]

    run_cached(sample_fraction = 0.5, aligner_args = `--beta`)
    @test calls == ["full", "sequence_subset_001", "sequence_subset_001",
        "full", "sequence_subset_001", "full", "sequence_subset_001"]

    changed_calls = String[]
    changed_aligner = (input_fasta, output_fasta; logs_dir = nothing,
        run_label = "run", aligner_args = Cmd(String[])) -> begin
        push!(changed_calls, String(run_label))
        return copy_aligner(input_fasta, output_fasta;
            logs_dir, run_label, aligner_args)
    end
    run_cached(aligner_fn = changed_aligner, sample_fraction = 0.5)
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
        return Iduna.EPLI.RawScore(raw)
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
    @test result.epli.value == 50.0
    scores = Iduna.EPLI.DataFrame(Iduna.EPLI.CSV.File(result.scores_path))
    @test scores.normalization_score == fill(20.0, 2)
end
