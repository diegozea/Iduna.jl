using .ExampleFixtures

@testset "README examples" begin
    @testset "Julia quick start" begin
        mktempdir() do tmp
            kwargs = merge((mmseqs_db = "/path/to/mmseqs/db",),
                iduna_stage_kwargs())
            result = cd(tmp) do
                Iduna.iduna("P20963"; kwargs...)
            end
            expanded = Iduna.load_expanded_msa(result)

            @test result.input_id == "P20963"
            @test Iduna.ResultsValidation.nsequences(expanded) == 1
        end
    end

    @testset "Command line quick start" begin
        command = [
            "julia",
            "--project=.",
            "-m",
            "Iduna",
            "P20963",
            "--mmseqs-db",
            "/path/to/mmseqs/db",
        ]
        kwargs = Iduna._parse_app_args(iduna_args_after_module(command))

        @test kwargs[:id] == "P20963"
        @test kwargs[:mmseqs_db] == "/path/to/mmseqs/db"
        @test !haskey(kwargs, :workdir)
    end

    @testset "ThorAxe-only example" begin
        mktempdir() do tmp
            kwargs = merge((
                    no_expansion = true,
                    workdir = joinpath(tmp, "ENST00000362089_thoraxe")),
                iduna_stage_kwargs(; fail_on_expansion = true))
            thoraxe_only = Iduna.iduna("ENST00000362089.10"; kwargs...)
            seed = Iduna.load_seed_msa(thoraxe_only)

            @test isempty(thoraxe_only.expansions)
            @test Iduna.ResultsValidation.nsequences(seed) == 1
        end
    end
end
