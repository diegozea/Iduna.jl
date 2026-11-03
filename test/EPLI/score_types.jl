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
    @test E._score_result(score) === score
    @test E.comparable_positions_normalization(score).normalized_score == 75
    @test E.comparable_positions_normalization(counts).normalized_score == 75
    @test !E.comparable_positions_normalization(score).no_comparable_positions
    @test_throws ArgumentError E.comparable_positions_normalization(E.RawScore(1))
    zero = E.comparable_positions_normalization(E.IdentityCounts(0, 0))
    @test zero isa E.NormalizationResult
    @test zero.normalized_score === 0.0
    @test zero.no_comparable_positions
    @test !E.comparable_positions_normalization(E.IdentityCounts(0, 4)).no_comparable_positions
    @test_throws ArgumentError E.NormalizationResult(1, 0; no_comparable_positions = true)
    @test_throws ArgumentError E.NormalizationResult(0; no_comparable_positions = true)
    @test E._normalization_result(zero) === zero
    @test E._normalization_result((; normalized_score = 0, normalization_score = 0,
        no_comparable_positions = true)).no_comparable_positions
    @test E.raw_score(E.RawScore(-2)) === -2.0
    @test E.no_normalization(E.RawScore(-2)).normalized_score === -2.0
    @test E.self_reference_normalization(E.RawScore(15);
        reference_score = E.RawScore(10)).normalized_score === 150.0
    @test E.self_reference_normalization(E.RawScore(15);
        reference_score = E.RawScore(-10)).normalized_score === -150.0
    @test !E.requires_reference_score(E.no_normalization)
    @test E.requires_reference_score(E.self_reference_normalization)
    @test E._score_result((; raw_score = 1, matched_positions = missing,
        comparable_positions = missing)).counts === nothing
    @test_throws ArgumentError E._score_result((; raw_score = 1,
        matched_positions = 3, comparable_positions = 2))
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
    @test E._parse_hhsuite_query_segment("Q A 1 A 1").cols == 7:7
    positions, codes = E._get_codes(block)
    @test E._identity_counts_from_codes(positions, codes) == E.IdentityCounts(2, 2)
    @test E._get_codes(replace(block, "Confidence\n" => "")) == (positions, codes)
    @test E._hhsuite_query_indices("A.~", 1, 2) == [1, 0, 2]
    # Use the first displayed query sequence for position counts, as before; the
    # comparison row is after the consensus row, whose insert columns may differ.
    consensus = "Probab=99\nQ query 1 Aac 3\nQ Conse 1 A.. 1\n          |  \nT query 1 ABC 3\nConfidence\n"
    @test E._identity_counts_from_codes(E._get_codes(consensus)...) ==
          E.IdentityCounts(1, 3)
    @test E._identity_counts_from_codes(E._get_codes(block * "\n" * block)...) ==
          E.IdentityCounts(2, 2)
    no_matches = replace(block, "| |" => "   ")
    @test E._identity_counts_from_codes(E._get_codes(no_matches)...).comparable_positions ==
          2
    @test E._identity_counts_from_codes(E._get_codes(no_matches)...).matched_positions == 0
    no_positions = replace(block, "A-C 2" => "--- 0", "| |" => "   ")
    @test E._identity_counts_from_codes(E._get_codes(no_positions)...) ==
          E.IdentityCounts(0, 0)
    @test E.comparable_positions_normalization(
        E._identity_counts_from_codes(E._get_codes(no_positions)...)).no_comparable_positions
    for malformed in ("", "unrecognized output", "Probab=99\n",
        replace(block, "Q query" => "query"), replace(block, "T query" => "target"),
        replace(block, "          | |\n" => ""),
        replace(block, "ABC 3" => "AB 2"), replace(block, "          | |" => " |"),
        replace(block, "| |" => "???"), block * "Q query 3 AA 4\n",
        replace(block, "Confidence\n" => "") * block)
        @test_throws ArgumentError E._get_codes(malformed)
    end
    @test_throws ArgumentError E._identity_counts_from_codes([1], Char[])
    @test_throws ArgumentError E._identity_counts_from_codes([-1], ['|'])
end

@testset "Optional EPLI and absence metadata" begin
    score = Iduna.EPLIScore(-2)
    @test Iduna.require_epli(score) === -2.0
    @test_throws ArgumentError Iduna.require_epli(nothing)
    seed_kwargs = (; pid = 10, stockholm_path = "seed.sto", summary_path = "summary.csv")
    for (epli, reason) in ((score, nothing), (nothing, nothing),
        (nothing, :sampling_disabled), (nothing, :not_recorded), (nothing, :ineligible))
        seed = Iduna.SeedSelection(; seed_kwargs..., epli, epli_reason = reason)
        @test seed.epli === epli
        @test seed.epli_reason === reason
        @test !hasproperty(seed, :score)
        relative = Iduna.Utils._relative_seed_paths(seed, pwd())
        @test relative.epli === epli
        @test relative.epli_reason === reason
        summary = Iduna.Utils._seed_summary(seed)
        @test !hasproperty(summary, :epli_state)
        data = Iduna.Utils.JSON.parse(Iduna.Utils.JSON.json(summary))
        loaded = Iduna._load_result_seed_record(data, pwd())
        @test loaded.epli === epli
        @test loaded.epli_reason === reason
    end
    @test Iduna.SeedSelection(; seed_kwargs..., epli = missing).epli === nothing
    @test Iduna.SeedSelection(; seed_kwargs..., epli = 150).epli.value === 150.0
    @test_throws ArgumentError Iduna.SeedSelection(; seed_kwargs...,
        epli = nothing, epli_reason = :unknown)
    @test_throws ArgumentError Iduna.SeedSelection(; seed_kwargs...,
        epli = score, epli_reason = :sampling_disabled)
    legacy = Dict("pid" => 10, "stockholm_path" => "seed.sto", "summary_path" => "s.csv")
    unavailable = Iduna._load_result_seed_record(legacy, pwd())
    @test unavailable.epli === nothing
    @test unavailable.epli_reason === :not_recorded
    @test Iduna._load_epli(Dict("epli" => 150)).value == 150
    @test Iduna._load_epli(Dict("epli" => nothing, "epli_state" => "not_computed")) ===
          nothing
    for data in (Dict("epli" => nothing, "epli_state" => "computed"),
        Dict("epli" => 2, "epli_state" => "not_computed"),
        Dict("epli" => NaN), Dict("epli" => 2, "epli_state" => "unknown"),
        Dict("epli" => 2, "epli_reason" => "sampling_disabled"),
        Dict("epli" => nothing, "epli_reason" => "unknown"))
        @test_throws ArgumentError Iduna._load_result_seed_record(merge(legacy, data), pwd())
    end
end

@testset "Sampled selection requires computed evidence" begin
    T = Iduna.ThorAxeMSA
    mktempdir() do tmp
        path = joinpath(tmp, "candidates.csv")
        header = "pid,pid_order,eligible,selected,epli,stockholm_path,fasta_path,selection_mode,pid_sample_count,epli_state,epli_reason\n"
        valid = "10,1,true,true,75,s10.sto,s10.fa,sampled_selection,2,computed,\n"
        for bad in (
            "20,2,true,false,,s20.sto,s20.fa,sampled_selection,2,not_computed,sampling_disabled\n",
            "20,2,true,false,,s20.sto,s20.fa,all_candidates,2,not_computed,sampling_disabled\n",
            "20,2,true,false,1,s20.sto,s20.fa,unknown,2,computed,\n",
            "20,2,true,false,1,s20.sto,s20.fa,sampled_selection,-1,computed,\n",
            "20,2,true,false,1,s20.sto,s20.fa,sampled_selection,2,not_computed,\n",
            "20,2,true,false,1,s20.sto,s20.fa,sampled_selection,2,computed,sampling_disabled\n",
            "20,2,true,false,,s20.sto,s20.fa,all_candidates,0,not_computed,ineligible\n")
            write(path, header * valid * bad)
            @test_throws ArgumentError T.select_best_seed(path)
            @test_throws ArgumentError T._selected_candidate_seeds(path)
        end
        write(path, header *
                    "10,1,true,true,,s10.sto,s10.fa,all_candidates,0,not_computed,sampling_disabled\n")
        @test only(T._selected_candidate_seeds(path)).epli_reason == :sampling_disabled
        @test_throws ArgumentError T.select_best_seed(path)
        write(path, header * replace(valid, "true,true" => "false,true"))
        @test_throws ArgumentError T._selected_candidate_seeds(path)
        write(path,
            "pid,pid_order,eligible,selected,epli,stockholm_path,fasta_path,selection_mode,pid_sample_count,epli_reason\n" *
            "10,1,true,true,,s10.sto,s10.fa,all_candidates,0,sampling_disabled\n")
        absent = only(T._selected_candidate_seeds(path))
        @test absent.epli === nothing
        @test absent.epli_reason === :sampling_disabled
        @test_throws ArgumentError T.select_best_seed(path)
        write(path, replace(read(path, String), "all_candidates,0" => "sampled_selection,2"))
        @test_throws ArgumentError T._selected_candidate_seeds(path)
        # Zero remains eligible; ties still use MSA size and then PID order.
        write(path,
            "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
            "10,2,true,0,2,s10.sto,s10.fa\n20,1,true,0,3,s20.sto,s20.fa\n")
        @test T.select_best_seed(path).pid == 20
        @test T._sampling_disabled((; pid_sample_count = 0))
        @test !T._sampling_disabled((; selection_mode = "sampled_selection"))
    end
end

@testset "Typed callbacks and zero-comparison diagnostics" begin
    E = Iduna.EPLI
    mktempdir() do tmp
        input = joinpath(tmp, "input.fa")
        write(input, ">ref\nACDE\n>other\nACDF\n")
        result = E.epli_score(input, joinpath(tmp, "typed"), copy_aligner;
            score_fn = TypedTestScore(), normalization_fn = TypedTestNormalizer(),
            sample_count = 2, progress_enabled = false)
        @test result.epli isa Iduna.EPLIScore
        @test result.epli.value == -20
        @test E.raw_score(result.reference_score) == 10
        @test result.no_comparable_samples == 0
        result = E.epli_score(input, joinpath(tmp, "zero"), copy_aligner;
            score_fn = (args...; kwargs...) -> E.IdentityScore(E.IdentityCounts(0, 0)),
            sample_count = 2, progress_enabled = false)
        @test result.epli.value === 0.0
        @test result.no_comparable_samples == 2
        @test all(row -> row.no_comparable_positions, E.CSV.File(result.scores_path))
        @test only(E.CSV.File(result.summary_path)).no_comparable_samples == 2
        T = Iduna.ThorAxeMSA
        @test_throws ArgumentError T._candidate_summary_row(
            Iduna.ResolvedTarget(input_id = "x", input_kind = :uniprot,
                ensembl_gene_id = "g", transcript_id = "t"), nothing, 10, 1,
            (; eligible = true); sample_count = 2, sample_fraction = 1,
            sample_seed = UInt64(1), metadata = nothing)
    end
end
