@testset "validation relative workdirs and relocation" begin
    mktempdir() do tmp
        target = Iduna.ResolvedTarget(; input_id = "ENST000001",
            input_kind = :ensembl_transcript, ensembl_gene_id = "GENE",
            transcript_id = "ENST000001")
        cd(tmp) do
            for path_kind in (:relative, :absolute, :dot)
                @testset "$path_kind" begin
                    root = joinpath(tmp, String(path_kind))
                    mkpath(root)
                    write(joinpath(root, "seed.sto"),
                        "# STOCKHOLM 1.0\nGENE ACDE\nORTHO ACDF\n//\n")
                    seed = Iduna.SeedSelection(; pid = 10.0, epli = 100.0,
                        stockholm_path = "seed.sto", summary_path = "summary.csv")
                    cd(path_kind === :dot ? root : tmp) do
                        requested = path_kind === :relative ? String(path_kind) :
                                    path_kind === :dot ? "." : root
                        result = Iduna.ResultsValidation.validate_results(
                            target, seed, nothing, requested)
                        @test result.stats_path ==
                              joinpath(root, "validation", "pid_10.00", "stats.csv")
                        @test result.seed_nseq == 2
                        reused = Iduna.ResultsValidation.validate_results(
                            target, seed, nothing, requested)
                        @test reused.seed_nseq == 2
                        state = Iduna.Utils._read_stage_state(dirname(reused.stats_path))
                        @test state["action"] == "reuse"
                        @test Set(keys(state["output_digests"])) ==
                              Set([joinpath("validation", "pid_10.00", "stats.csv")])
                    end
                    relocated = joinpath(tmp, "moved_$(path_kind)")
                    mv(root, relocated)
                    result = Iduna.ResultsValidation.validate_results(
                        target, seed, nothing, basename(relocated))
                    @test result.seed_nseq == 2
                    @test startswith(result.stats_path, relocated)
                    @test Iduna.Utils._read_stage_state(dirname(result.stats_path))["action"] ==
                          "reuse"
                end
            end
        end
    end
end

@testset "nested alias cannot certify corrupted validation output" begin
    mktempdir() do tmp
        cd(tmp) do
            mkpath("out")
            write(joinpath("out", "seed.sto"),
                "# STOCKHOLM 1.0\nGENE ACDE\nORTHO ACDF\n//\n")
            target = Iduna.ResolvedTarget(; input_id = "ENST000001",
                input_kind = :ensembl_transcript, ensembl_gene_id = "GENE",
                transcript_id = "ENST000001")
            seed = Iduna.SeedSelection(; pid = 10.0, epli = 100.0,
                stockholm_path = "seed.sto", summary_path = "summary.csv")
            result = Iduna.ResultsValidation.validate_results(target, seed, nothing, abspath("out"))
            actual = joinpath("out", "validation", "pid_10.00", "stats.csv")
            alias = joinpath("out", actual)
            mkpath(dirname(alias))
            cp(actual, alias)
            # Give the relative call a valid nested copy before corrupting the real file.
            Iduna.ResultsValidation.validate_results(target, seed, nothing, "out")
            table = Iduna.ResultsValidation.DataFrame(Iduna.ResultsValidation.CSV.File(actual))
            table.seed_nseq .= 99
            Iduna.ResultsValidation.CSV.write(actual, table)
            corrected = Iduna.ResultsValidation.validate_results(target, seed, nothing, "out")
            @test corrected.seed_nseq == 2
            @test Iduna.Utils._read_stage_state(dirname(actual))["action"] == "rebuild"
            @test read(alias, String) == read(corrected.stats_path, String)
        end
    end
end
