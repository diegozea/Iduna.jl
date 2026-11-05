struct TypedTestScore end
function (::TypedTestScore)(reference, sample; label, kwargs...)
    Iduna.EPLI.RawScore(label == "reference_self" ? 10 : -2)
end
struct TypedTestNormalizer end
function (::TypedTestNormalizer)(score; reference_score)
    Iduna.EPLI.self_reference_normalization(score; reference_score)
end
Iduna.EPLI.requires_reference_score(::TypedTestNormalizer) = true

@testset "Scientific score invariants" begin
    E = Iduna.EPLI
    for invalid in (NaN, Inf, -Inf, big"1e1000")
        @test_throws ArgumentError E.RawScore(invalid)
        @test_throws ArgumentError E.NormalizationResult(invalid)
        @test_throws ArgumentError E.NormalizationResult(1, invalid)
        @test_throws ArgumentError Iduna.EPLIScore(invalid)
        @test_throws ArgumentError Iduna.SeedSelection(pid = 10, epli = invalid,
            stockholm_path = "seed.sto", summary_path = "summary.csv")
    end
    for pair in ((-1, 2), (3, 2), (1, 0), (0.5, 2), (0, Inf), (NaN, 2),
        (0, big(typemax(Int)) + 1))
        @test_throws ArgumentError E.IdentityCounts(pair...)
    end
    counts = E.IdentityCounts(3.0, 4)
    score = E.IdentityScore(counts; diagnostics = (; probability = missing))
    @test E.raw_score(score) === 3.0
    @test score isa E.ScoreResult
    @test E.comparable_positions_normalization(score).normalized_score == 75
    @test_throws ArgumentError E.comparable_positions_normalization(E.RawScore(1))
    zero = E.comparable_positions_normalization(E.IdentityScore(E.IdentityCounts(0, 0)))
    @test zero isa E.NormalizationResult
    @test zero.normalized_score === 0.0
    @test zero.normalization_score === 0.0
    @test E.comparable_positions_normalization(E.IdentityScore(E.IdentityCounts(0, 4))).normalization_score ==
          4
    for value in (0.0, -0.0), denominator in (0.0, -0.0)

        result = E.NormalizationResult(value, denominator)
        @test result.normalized_score === value
        @test result.normalization_score === denominator
    end
    @test E.raw_score(E.RawScore(-2)) === -2.0
    @test E.no_normalization(E.RawScore(-2)).normalized_score === -2.0
    @test E.self_reference_normalization(E.RawScore(15);
        reference_score = E.RawScore(10)).normalized_score === 150.0
    @test E.self_reference_normalization(E.RawScore(15);
        reference_score = E.RawScore(-10)).normalized_score === -150.0
    @test !E.requires_reference_score(E.no_normalization)
    @test E.requires_reference_score(E.self_reference_normalization)
    @test E.RawScore(1) isa E.ScoreResult
    @test E.RawScore(1).counts === nothing
    diagnostics = Dict(:value => NaN)
    raw = E.RawScore(2; diagnostics)
    diagnostics[:value] = Inf
    @test E.raw_score(raw) == 2 # Mutable diagnostics cannot alter the score invariant.
    @test !ismutable(raw)
    @test !ismutable(counts)
end

@testset "Recognized HHsuite evidence" begin
    E = Iduna.EPLI
    block = "Probab=99\nQ query 1 A-C 2\n          | |\nT query 1 ABC 3\nConfidence\n"
    @test E._parse_hhsuite_sequence_row("Q A 1 A 1").cols == 7:7
    positions, codes = E._get_codes(block)
    @test E._identity_counts_from_codes(positions, codes) == E.IdentityCounts(2, 2)
    @test E._get_codes(replace(block, "Confidence\n" => "")) == (positions, codes)
    @test E._hhsuite_residue_indices("A.~", 1, 2) == [1, 0, 2]
    @test E._hhsuite_residue_indices("--.", 1, 0) == [0, 0, 0]
    @test_throws ArgumentError E._hhsuite_residue_indices("A", 0, 0)
    @test_throws ArgumentError E._hhsuite_residue_indices("A", -1, -1)
    @test_throws ArgumentError E._hhsuite_residue_indices("AA", 1, 1)
    annotated = replace(block, "Q query" => "Q ss_pred   HHH\nQ query")
    @test E._get_codes(annotated) == (positions, codes)
    # HHsuite accepts annotation names with these prefixes, including suffixes.
    for annotation in ("ss_pred", "ss_conf", "ss_dssp", "sa_dssp", "ss_pred_extra")
        symbols = annotation in ("ss_conf", "sa_dssp") ? "999" : "HHH"
        annotated = replace(block, "Q query" => "Q $annotation   $symbols\nQ query",
            "Confidence" => "T $annotation   $symbols\nConfidence")
        @test E._get_codes(annotated) == (positions, codes)
    end
    # Use the first displayed query sequence for position counts, as before; the
    # comparison row is after the consensus row, whose insert columns may differ.
    consensus = "Probab=99\nQ query 1 Aac 3\nQ Conse 1 A.. 1\n          |  \nT query 1 ABC 3\nConfidence\n"
    @test E._identity_counts_from_codes(E._get_codes(consensus)...) ==
          E.IdentityCounts(1, 3)
    # Validate all displayed sequences. A truncated first row must not fall back
    # to consensus counts; a damaged later row must not be silently ignored.
    for malformed in (
        replace(consensus, "Q query 1 Aac 3" => "Q query 1 Aac"),
        replace(consensus, "Q Conse 1 A.. 1" => "Q Conse 1 A.."),
        replace(consensus, "Q Conse 1 A.. 1" => "Q Conse 1 A. 1"),
        replace(consensus, "T query 1 ABC 3" => "T query 1 ABC\nT Conse 1 ABC 3"),
        replace(consensus, "T query 1 ABC 3" => "T query 1 ABC 3\nT other 1 ABC"),
        replace(consensus, "T query 1 ABC 3" => "T query 1 ABC 3\nT other 1 AB 2"),
        replace(consensus, "Q query 1 Aac 3" => "Q ss_pred 1 Aac"))
        @test_throws ArgumentError E._get_codes(malformed)
    end
    @test E._identity_counts_from_codes(E._get_codes(block * "\n" * block)...) ==
          E.IdentityCounts(2, 2)
    no_matches = replace(block, "| |" => "   ")
    @test E._identity_counts_from_codes(E._get_codes(no_matches)...).comparable_positions ==
          2
    @test E._identity_counts_from_codes(E._get_codes(no_matches)...).matched_positions == 0
    no_positions = replace(block, "A-C 2" => "--- 0", "| |" => "   ")
    @test E._identity_counts_from_codes(E._get_codes(no_positions)...) ==
          E.IdentityCounts(0, 0)
    @test E.comparable_positions_normalization(E.IdentityScore(
        E._identity_counts_from_codes(E._get_codes(no_positions)...))).normalized_score ===
          0.0
    for malformed in ("", "unrecognized output", "Probab=99\n",
        replace(block, "Q query" => "query"), replace(block, "T query" => "target"),
        replace(block, "          | |\n" => ""),
        replace(block, "ABC 3" => "AB 2"), replace(block, "          | |" => " |"),
        replace(block, "| |" => "???"), block * "Q query 3 AA 4\n",
        replace(block, "Confidence\n" => "") * block,
        replace(block, "Q query 1 A-C 2" => "Q query 0 A-C 1"),
        replace(block, "T query 1 ABC 3" => "T query 0 ABC 2"),
        replace(block, "Q query 1 A-C 2" => "Q query -1 A-C 0"),
        replace(block, "Q query 1 A-C 2" => "Q ss_pred   HHH"),
        replace(block, "T query 1 ABC 3" => "T ss_pred   HHH"))
        @test_throws ArgumentError E._get_codes(malformed)
    end
    @test_throws ArgumentError E._identity_counts_from_codes([1], Char[])
    @test_throws ArgumentError E._identity_counts_from_codes([-1], ['|'])
end

@testset "Complete HHsuite hit sections" begin
    E = Iduna.EPLI
    # PrintHHR uses No/name/Probab headers for each alternative hit; a hit can
    # contain several alignment blocks. Keep the native sequence-column layout.
    output = """
    No 1
    >query
    Probab=99.00  E-value=1e-10  Score=20.00

    Q query             1 AC    2 (4)
    Q Consensus         1 ac    2 (4)
                          |+
    T Consensus         1 ac    2 (4)
    T query             1 AC    2 (4)
    Confidence            99


    Q query             3 DE    4 (4)
    Q Consensus         3 de    4 (4)
                          .|
    T Consensus         3 de    4 (4)
    T query             3 DE    4 (4)
    Confidence            99


    No 2
    >query
    Probab=98.00  E-value=1e-9  Score=18.00

    Q query             1 AC    2 (4)
    Q Consensus         1 ac    2 (4)
                          ||
    T Consensus         3 ac    4 (4)
    T query             3 AC    4 (4)
    Confidence            99

    """
    positions, codes = E._get_codes(output)
    @test positions == [1, 2, 3, 4, 1, 2]
    @test E._identity_counts_from_codes(positions, codes) == E.IdentityCounts(3, 4)
    @test E._get_codes(replace(output, "Confidence            99\n" => "")) ==
          (positions, codes)
    first_hit, second_hit = split(output, "No 2\n"; limit = 2)
    second_header = first(split(second_hit, "Q query"; limit = 2))
    for incomplete in (
        "Probab=99\n\n" * output,
        first_hit * "Probab=99\n\nNo 2\n" * second_hit,
        first_hit * "No 2\n" * second_header,
        output * "No 3\n>query\nProbab=99\n\n",
        "No 0\n>query\n" * output,
        first_hit * "No 0\n>query\nNo 2\n" * second_hit,
        first_hit * "No 2\n",
        first_hit * "No 2\n>query\n",
        first_hit * "No 2\n>query\nProbab",
        replace(output, "Probab=99.00  E-value=1e-10  Score=20.00\n" => ""),
        replace(output, "Probab=98.00  E-value=1e-9  Score=18.00\n" => ""))
        @test_throws ArgumentError E._get_codes(incomplete)
    end
end

@testset "HHsuite numeric fields" begin
    E = Iduna.EPLI
    summary = E._hhalign_summary("Probab=99.00  E-value=1.2e-30  Score=-3.50")
    @test summary == (; probability = 99.0, e_value = 1.2e-30, hhalign_score = -3.5)
    @test E._parse_hhalign_float("OtherScore=10 Score=2", "Score") == 2.0
    @test ismissing(E._parse_hhalign_float("OtherScore=10", "Score"))
    for token in ("1e+", "1.2garbage", "1e999", "")
        @test ismissing(E._parse_hhalign_float("Score=$(token)", "Score"))
        diagnostics = E._hhalign_summary("Probab=99 E-value=$(token) Score=10")
        @test ismissing(diagnostics.e_value)
        @test E.raw_score(E.IdentityScore(E.IdentityCounts(1, 2); diagnostics)) == 1.0
        @test E.raw_score(E.RawScore(diagnostics.hhalign_score; diagnostics)) == 10.0
    end
    for (token, value) in (("nan", NaN), ("-nan", NaN), ("inf", Inf),
        ("+inf", Inf), ("-inf", -Inf))
        diagnostics = E._hhalign_summary("Probab=99 E-value=$(token) Score=$(token)")
        @test isequal(diagnostics.e_value, value)
        @test isequal(diagnostics.hhalign_score, value)
        @test E.raw_score(E.IdentityScore(E.IdentityCounts(1, 2); diagnostics)) == 1.0
        @test_throws ArgumentError E.RawScore(diagnostics.hhalign_score)
    end
end

@testset "Optional EPLI and absence metadata" begin
    score = Iduna.EPLIScore(-2)
    seed_kwargs = (; pid = 10, stockholm_path = "seed.sto", summary_path = "summary.csv")
    for (epli, reason) in ((score, nothing), (nothing, nothing),
        (nothing, :sampling_disabled), (nothing, :not_recorded), (nothing, :ineligible))
        seed = Iduna.SeedSelection(; seed_kwargs..., epli, epli_reason = reason)
        @test seed.epli === epli
        @test seed.epli_reason === reason
        if epli === nothing
            @test_throws ArgumentError Iduna.epli_value(seed)
        else
            @test Iduna.epli_value(seed) === -2.0
        end
        @test !hasproperty(seed, :score)
        relative = Iduna.Utils._relative_seed_paths(seed, pwd())
        @test relative.epli === epli
        @test relative.epli_reason === reason
        summary = Iduna.Utils._seed_summary(seed)
        data = Iduna.Utils.JSON.parse(Iduna.Utils.JSON.json(summary))
        loaded = Iduna._load_result_seed_record(data, pwd())
        @test loaded.epli === epli
        @test loaded.epli_reason === reason
    end
    @test Iduna.SeedSelection(; seed_kwargs..., epli = missing).epli === nothing
    @test Iduna.epli_value(Iduna.SeedSelection(; seed_kwargs..., epli = 150)) === 150.0
    @test_throws ArgumentError Iduna.SeedSelection(; seed_kwargs...,
        epli = nothing, epli_reason = :unknown)
    @test_throws ArgumentError Iduna.SeedSelection(; seed_kwargs...,
        epli = score, epli_reason = :sampling_disabled)
    legacy = Dict("pid" => 10, "stockholm_path" => "seed.sto", "summary_path" => "s.csv")
    unavailable = Iduna._load_result_seed_record(legacy, pwd())
    @test unavailable.epli === nothing
    @test unavailable.epli_reason === :not_recorded
    @test Iduna._load_result_seed_record(merge(legacy, Dict("epli" => 150)), pwd()).epli.value ==
          150
    for data in (Dict("epli" => NaN), Dict("epli" => Inf),
        Dict("epli" => 2, "epli_reason" => "sampling_disabled"),
        Dict("epli" => nothing, "epli_reason" => "unknown"))
        @test_throws ArgumentError Iduna._load_result_seed_record(merge(legacy, data), pwd())
    end
end

@testset "Sampled selection requires scores" begin
    T = Iduna.ThorAxeMSA
    mktempdir() do tmp
        path = joinpath(tmp, "candidates.csv")
        header = "pid,pid_order,eligible,selected,epli,stockholm_path,fasta_path,selection_mode,pid_sample_count,epli_reason\n"
        valid = "10,1,true,true,75,s10.sto,s10.fa,sampled_selection,2,\n"
        for bad in (
            "20,2,true,false,,s20.sto,s20.fa,sampled_selection,2,sampling_disabled\n",
            "20,2,true,false,,s20.sto,s20.fa,all_candidates,2,sampling_disabled\n",
            "20,2,true,false,1,s20.sto,s20.fa,unknown,2,\n",
            "20,2,true,false,1,s20.sto,s20.fa,sampled_selection,-1,\n",
            "20,2,true,false,invalid,s20.sto,s20.fa,sampled_selection,2,\n",
            "20,2,true,false,1,s20.sto,s20.fa,sampled_selection,2,sampling_disabled\n",
            "20,2,true,false,,s20.sto,s20.fa,all_candidates,0,ineligible\n")
            write(path, header * valid * bad)
            @test_throws ArgumentError T.select_best_seed(path)
            @test_throws ArgumentError T._selected_candidate_seeds(path)
        end
        # Absence reasons are optional when sampling metadata establishes absence.
        for reason in ("sampling_disabled", "")
            write(path, header *
                        "10,1,true,true,,s10.sto,s10.fa,all_candidates,0,$reason\n")
            absent = only(T._selected_candidate_seeds(path))
            @test absent.epli === nothing
            @test absent.epli_reason === :sampling_disabled
            @test_throws ArgumentError T.select_best_seed(path)
        end
        write(path, header * replace(valid, "true,true" => "false,true"))
        @test_throws ArgumentError T._selected_candidate_seeds(path)
        write(path, header * replace(valid, "true,true" => "true,false"))
        @test isempty(T._selected_candidate_seeds(path))

        # Filtering and ordering must keep each parsed score attached to its PID.
        write(path,
            header *
            "30,3,true,false,-2,s30.sto,s30.fa,sampled_selection,2,\n" * valid *
            "20,2,false,false,NaN,s20.sto,s20.fa,sampled_selection,2,\n" *
            "40,4,true,true,150,s40.sto,s40.fa,sampled_selection,2,\n")
        best = T.select_best_seed(path)
        @test (best.pid, best.epli.value) == (40, 150)
        @test [(s.pid, s.epli.value) for s in T._selected_candidate_seeds(path)] ==
              [(10, 75), (40, 150)]
        @test [(s.pid, s.epli.value) for s in T._eligible_candidate_seeds(path)] ==
              [(10, 75), (30, -2), (40, 150)]
        # Without eligibility metadata, every candidate must be valid before cache reuse.
        write(path,
            "pid,selected,epli,stockholm_path,fasta_path\n" *
            "10,true,75,s10.sto,s10.fa\n20,false,NaN,s20.sto,s20.fa\n")
        @test_throws ArgumentError T._selected_candidate_seeds(path)
        # Zero remains eligible; ties still use MSA size and then PID order.
        write(path,
            "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
            "10,2,true,0,2,s10.sto,s10.fa\n20,1,true,0,3,s20.sto,s20.fa\n")
        @test T.select_best_seed(path).pid == 20
        # Signed zeros are the same scientific score. Both MSA size and PID order
        # must break their ties, while the selected seed retains its original sign.
        for (first_score, second_score) in (("-0.0", "0.0"), ("0.0", "-0.0"))
            for (first_size, second_size) in ((5, 2), (3, 3))
                write(path,
                    "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
                    "10,1,true,$first_score,$first_size,s10.sto,s10.fa\n" *
                    "20,2,true,$second_score,$second_size,s20.sto,s20.fa\n")
                best = T.select_best_seed(path)
                @test best.pid == 10
                @test best.epli.value === parse(Float64, first_score)
            end
        end
        # Rank the validated Float64 values. These integer inputs become the same
        # score, so the larger MSA must win the tie.
        write(path,
            "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
            "10,2,true,9007199254740992,5,s10.sto,s10.fa\n" *
            "20,1,true,9007199254740993,2,s20.sto,s20.fa\n")
        best = T.select_best_seed(path)
        @test best.pid == 10
        @test best.epli.value == 9007199254740992.0
        @test T._sampling_disabled((; pid_sample_count = 0))
        @test !T._sampling_disabled((; selection_mode = "sampled_selection"))
        @test T._sampling_disabled((; selection_mode = "all_candidates"))
    end
end

@testset "Typed callbacks and zero-comparison diagnostics" begin
    E = Iduna.EPLI
    mktempdir() do tmp
        input = write_epli_test_input(tmp)
        result = E.epli_score(input, joinpath(tmp, "typed"), copy_aligner;
            score_fn = TypedTestScore(), normalization_fn = TypedTestNormalizer(),
            sample_count = 2, progress_enabled = false)
        @test result.epli isa Iduna.EPLIScore
        @test result.epli.value == -20
        @test E.raw_score(result.reference_score) == 10
        @test result.no_comparable_samples == 0
        @test result.samples_without_identity_counts == 2
        @test all(row -> ismissing(row.no_comparable_positions), result.rows)
        for normalization_fn in (E.comparable_positions_normalization,
            E.no_normalization, E.self_reference_normalization,
            (score; kwargs...) -> E.NormalizationResult(0))
            score_fn = (args...; label, kwargs...) -> E.IdentityScore(
                label == "reference_self" ? E.IdentityCounts(2, 2) : E.IdentityCounts(0, 0))
            result = E.epli_score(
                input, joinpath(tmp, string(normalization_fn)), copy_aligner;
                score_fn, normalization_fn, sample_count = 2, progress_enabled = false)
            @test result.epli.value === 0.0
            @test result.no_comparable_samples == 2
            @test result.samples_without_identity_counts == 0
            @test all(row -> row.no_comparable_positions, E.CSV.File(result.scores_path))
            @test only(E.CSV.File(result.summary_path)).no_comparable_samples == 2
        end
        # The evidence follows validated counts even for negative profile scores,
        # and unavailable counts remain distinct from observed zero matches.
        scores = [E.RawScore(-2; counts = E.IdentityCounts(0, 0)),
            E.IdentityScore(E.IdentityCounts(0, 4)), E.RawScore(0)]
        score_fn = (args...; label, kwargs...) -> scores[parse(Int, last(split(label, '_')))]
        result = E.epli_score(input, joinpath(tmp, "mixed_counts"), copy_aligner;
            score_fn, normalization_fn = E.no_normalization,
            sample_count = 3, progress_enabled = false)
        @test isequal([row.no_comparable_positions for row in result.rows], [
            true, false, missing])
        @test result.epli.value == 0
        @test result.rows[1].normalized_score == -2
        @test result.no_comparable_samples == 1
        @test result.samples_without_identity_counts == 1
        @test only(E.CSV.File(result.summary_path)).samples_without_identity_counts == 1
        T = Iduna.ThorAxeMSA
        @test_throws ArgumentError T._candidate_summary_row(
            Iduna.ResolvedTarget(input_id = "x", input_kind = :uniprot,
                ensembl_gene_id = "g", transcript_id = "t"), nothing, 10, 1,
            (; eligible = true); sample_count = 2, sample_fraction = 1,
            sample_seed = UInt64(1), metadata = nothing)
    end
end
