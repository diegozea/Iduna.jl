@testset "Normalization at Float64 limits" begin
    normalize = (raw, reference) -> Iduna.EPLI.self_reference_normalization(
        Iduna.EPLI.RawScore(raw); reference_score = Iduna.EPLI.RawScore(reference)).normalized_score
    tiny = nextfloat(0.0)
    largest = floatmax(Float64)
    for (raw, reference, expected) in (
        (1e308, 1e308, 100.0),
        (-1e308, 1e308, -100.0),
        (1e308, -1e308, -100.0),
        (-1e308, -1e308, 100.0),
        (largest, 100.0, largest),
        (prevfloat(largest), 100.0, prevfloat(largest)),
        (9.047102574464534e307, 50.32617858418559, largest),
        (tiny, tiny, 100.0),
        (tiny, 100.0, tiny),
        (tiny, 1.0, 100 * tiny),
        (tiny, 200.0, 0.0),
        (3 * tiny, 200.0, 2 * tiny),
        (-tiny, 200.0, -0.0),
        (tiny, largest, 0.0),
        (0.0, -1.0, -0.0),
        (-0.0, -1.0, 0.0),
        (15.0, 10.0, 150.0))
        @test isequal(normalize(raw, reference), expected)
    end
    for raw in (largest, -largest)
        @test_throws ArgumentError normalize(raw, 1.0)
        @test_throws ArgumentError normalize(raw, tiny)
    end
    for reference in (0.0, -0.0)
        @test_throws ErrorException normalize(1.0, reference)
    end
    # Numerical extremes belong to generic ratios, not impossible identity counts.
    @test Iduna.EPLI._percentage_ratio(1e308, 100.0) == 1e308
    @test Iduna.EPLI._percentage_ratio(tiny, 100.0) == tiny
    @test Iduna.EPLI.comparable_positions_normalization(
        Iduna.EPLI.IdentityScore(Iduna.EPLI.IdentityCounts(typemax(Int), typemax(Int)))).normalized_score ==
          100.0

    # An exact rational oracle checks combinations across the Float64 exponent range.
    values = (tiny, 3 * tiny, floatmin(Float64), 1e-300, 0.125, 3.0, 1e300, largest)
    for raw in values, reference in values

        exact = 100 * Rational{BigInt}(raw) / Rational{BigInt}(reference)
        expected = Float64(BigFloat(exact, RoundNearest; precision = 512), RoundNearest)
        if isfinite(expected)
            @test isapprox(normalize(raw, reference), expected;
                atol = eps(expected), rtol = 2 * eps(Float64))
        else
            @test_throws ArgumentError normalize(raw, reference)
        end
    end

    tasks = [Threads.@spawn setprecision(BigFloat, bits) do
                 setrounding(BigFloat, RoundUp) do
                     before = (precision(BigFloat), rounding(BigFloat))
                     results = [normalize(1e308, 1e308), normalize(tiny, 100.0)]
                     return (; before, after = (precision(BigFloat), rounding(BigFloat)),
                         results)
                 end
             end for bits in (32, 128)]
    for task in tasks
        result = fetch(task)
        @test result.before == result.after
        @test result.results == [100.0, tiny]
    end
end

@testset "Subnormal medians" begin
    tiny = nextfloat(0.0)
    @test Iduna.EPLI._score_median([tiny, tiny]) == tiny
    @test Iduna.EPLI._score_median([tiny, 2 * tiny]) == 2 * tiny
    @test Iduna.EPLI._score_median([-tiny, -tiny]) == -tiny
    @test Iduna.EPLI._score_median([-tiny, tiny]) == 0.0
    @test Iduna.EPLI._score_median([0.0, tiny]) == 0.0
    @test Iduna.EPLI._score_median([tiny, tiny, tiny]) == tiny
    @test Iduna.EPLI._score_median([floatmax(Float64), floatmax(Float64)]) ==
          floatmax(Float64)
end

@testset "Extreme normalization through scoring callbacks" begin
    mktempdir() do tmp
        input = write_epli_test_input(tmp)
        for (raw, reference, expected) in ((1e308, 1e308, 100.0), (
            nextfloat(0.0), 100.0, nextfloat(0.0)))
            score_fn = (reference_msa, sample_msa; label, kwargs...) -> begin
                value = label == "reference_self" ? reference : raw
                return Iduna.EPLI.RawScore(value)
            end
            result = Iduna.EPLI.epli_score(input, joinpath(tmp, string(raw)), copy_aligner;
                score_fn, normalization_fn = Iduna.EPLI.self_reference_normalization,
                sample_count = 2, sample_seed = 1, progress_enabled = false)
            @test result.epli.value == expected
            @test all(row.normalized_score == expected for row in result.rows)
            @test isfile(result.scores_path)
            @test isfile(result.summary_path)
        end
    end
end
