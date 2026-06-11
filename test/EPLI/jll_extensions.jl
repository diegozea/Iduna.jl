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
