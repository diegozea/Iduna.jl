using .ExampleFixtures
using MIToS.MSA: Stockholm, getannotcolumn, getannotfile, read_file

@testset "Docs examples" begin
    function _docs_full_result(tmp::AbstractString)
        kwargs = merge((
                mmseqs_db = "/path/to/mmseqs/uniref_db",
                workdir = joinpath(tmp, "P20963"),
                overwrite = false,
                centroids = false),
            iduna_stage_kwargs())
        return Iduna.iduna("P20963"; kwargs...)
    end

    @testset "home page Julia API example" begin
        mktempdir() do tmp
            result = _docs_full_result(tmp)
            expanded = Iduna.load_expanded_msa(result)

            @test result.input_id == "P20963"
            @test result.status === :ok
            @test Iduna.ResultsValidation.nsequences(expanded) == 1
        end
    end

    @testset "home page command-line examples" begin
        source_command = [
            "julia",
            "--project=.",
            "-m",
            "Iduna",
            "P20963",
            "--mmseqs-db",
            "/path/to/mmseqs/uniref_db",
        ]
        source_kwargs = Iduna._parse_app_args(iduna_args_after_module(source_command))
        @test source_kwargs[:id] == "P20963"
        @test source_kwargs[:mmseqs_db] == "/path/to/mmseqs/uniref_db"

        thoraxe_only_command = [
            "julia",
            "--project=.",
            "-m",
            "Iduna",
            "ENST00000362089.10",
            "--no-expansion",
        ]
        thoraxe_kwargs = Iduna._parse_app_args(
            iduna_args_after_module(thoraxe_only_command))
        @test thoraxe_kwargs[:id] == "ENST00000362089.10"
        @test thoraxe_kwargs[:no_expansion] === true
        @test !haskey(thoraxe_kwargs, :mmseqs_db)
    end

    @testset "home page ThorAxe-only Julia example" begin
        mktempdir() do tmp
            kwargs = merge((
                    no_expansion = true,
                    workdir = joinpath(tmp, "ENST00000362089_thoraxe")),
                iduna_stage_kwargs(; fail_on_expansion = true))
            thoraxe_only = Iduna.iduna("ENST00000362089.10"; kwargs...)

            seed = Iduna.load_seed_msa(thoraxe_only)
            baseline_stockholm = thoraxe_only.thoraxe_msa.baseline_stockholms[1]
            seed_stockholm = thoraxe_only.thoraxe_msa.seeds[1].stockholm_path

            @test isempty(thoraxe_only.expansions)
            @test Iduna.ResultsValidation.nsequences(seed) == 1
            @test isfile(artifact_path(baseline_stockholm, thoraxe_only.workdir))
            @test isfile(artifact_path(seed_stockholm, thoraxe_only.workdir))
        end
    end

    @testset "home page threading command examples" begin
        env_app_command = [
            "JULIA_NUM_THREADS=4",
            "iduna",
            "P20963",
            "--mmseqs-db",
            "/path/to/mmseqs/uniref_db",
            "--mmseqs-threads",
            "8",
        ]
        env_kwargs = Iduna._parse_app_args(iduna_args_after_app(env_app_command))
        @test env_kwargs[:id] == "P20963"
        @test env_kwargs[:mmseqs_threads] == 8

        separator_command = [
            "iduna",
            "--threads=4",
            "--",
            "P20963",
            "--mmseqs-db",
            "/path/to/mmseqs/uniref_db",
            "--mmseqs-threads",
            "8",
        ]
        separator_kwargs = Iduna._parse_app_args(
            iduna_args_after_app(separator_command))
        @test separator_kwargs[:id] == "P20963"
        @test separator_kwargs[:mmseqs_threads] == 8

        threaded_source_command = [
            "julia",
            "--threads",
            "4",
            "--project=.",
            "-m",
            "Iduna",
            "P20963",
            "--mmseqs-db",
            "/path/to/mmseqs/uniref_db",
            "--mmseqs-threads",
            "8",
        ]
        source_kwargs = Iduna._parse_app_args(
            iduna_args_after_module(threaded_source_command))
        @test source_kwargs[:id] == "P20963"
        @test source_kwargs[:mmseqs_threads] == 8
    end

    @testset "output page result path examples" begin
        mktempdir() do tmp
            result = _docs_full_result(tmp)

            @test isfile(artifact_path(
                result.thoraxe_msa.baseline_stockholms[1], result.workdir))
            @test isfile(artifact_path(
                result.thoraxe_msa.baseline_fastas[1], result.workdir))
            @test isfile(artifact_path(
                result.thoraxe_msa.seeds[1].stockholm_path, result.workdir))
            @test isfile(artifact_path(
                result.thoraxe_msa.seeds[1].fasta_path, result.workdir))
            @test isfile(artifact_path(
                result.thoraxe_msa.seeds[1].s_exon_blocks_tsv, result.workdir))
            @test result.thoraxe_msa.pid_sample_count == 1
            @test result.thoraxe_msa.pid_sample_fraction == 1.0
            @test result.thoraxe_msa.pid_sample_seed == UInt64(7)
            @test result.thoraxe_msa.sampling_strategy === :common

            expansion = result.expansions[1]
            @test !ismissing(expansion)
            @test isfile(artifact_path(expansion.match_stockholm, result.workdir))
            @test isfile(artifact_path(expansion.s_exon_blocks_tsv, result.workdir))
            @test isfile(artifact_path(expansion.a3m_path, result.workdir))
        end
    end

    @testset "output page block table examples" begin
        mktempdir() do tmp
            result = _docs_full_result(tmp)

            seed_blocks = joinpath(
                result.workdir,
                result.thoraxe_msa.seeds[1].s_exon_blocks_tsv,
            )
            expanded_blocks = joinpath(
                result.workdir,
                something(result.expansions[1]).s_exon_blocks_tsv,
            )

            @test isfile(seed_blocks)
            @test isfile(expanded_blocks)
        end
    end

    @testset "output page MIToS annotation example" begin
        mktempdir() do tmp
            stockholm_path = annotated_stockholm(joinpath(tmp, "annotated.sto"))

            msa = read_file(stockholm_path, Stockholm; keepinserts = true)
            column_codes = getannotcolumn(msa, "SExonCode")
            code_key = getannotfile(msa, "SExonCodeMap")

            @test column_codes == "aabb"
            @test occursin("\"a\"=>\"1_0\"", code_key)
            @test occursin("\"b\"=>\"12_2\"", code_key)
        end
    end
end
