@testset "protein overrides retain only available sequence artifacts" begin
    function override_resolver(calls; fetch_override = false, damage = nothing)
        proteins = fetch_override ? ["ENSP000001", "ENSP000002"] : ["ENSP000001"]
        transcripts = ["ENST00000$(index)" for index in eachindex(proteins)]
        entry = Iduna.IDMapping._UniProtEntry(;
            id = "P20963", species = "Homo sapiens", gene_ids = ["ENSG000001"],
            transcript_ids = transcripts,
            transcript_to_protein = Dict(zip(transcripts, proteins)),
            transcript_to_gene = Dict(transcript => "ENSG000001"
            for transcript in transcripts),
            transcript_to_isoform = Dict{String, String}(),
            protein_sequence = "ACDE")
        return (primary; kwargs...) -> begin
            calls[] += 1
            Iduna.IDMapping.resolve_target(primary; kwargs...,
                _uniprot_entry_fetcher = _ -> entry,
                _candidate_validator = (fetched, sequence_dir) -> begin
                    candidates, uniprot_path = Iduna.IDMapping._validate_candidates(
                        fetched, sequence_dir; _ensembl_protein_fetcher = _ -> "ACDE")
                    damage === nothing || damage(sequence_dir)
                    return candidates, uniprot_path
                end)
        end
    end
    stage_options = (; supplied_uniprot = nothing, ensembl_gene_id = nothing,
        disambiguating_transcript = nothing, species = nothing, overwrite = false)

    @testset "relative resolver workdirs return readable artifact paths" begin
        mktempdir() do root
            cd(root) do
                calls = Ref(0)
                resolver = override_resolver(calls)
                target = resolver("P20963"; workdir = "relative_resolver")
                @test isabspath(target.workdir)
                @test realpath(target.workdir) ==
                      realpath(joinpath(root, "relative_resolver"))
                @test isabspath(target.uniprot_sequence_path)
                @test isabspath(target.ensembl_protein_sequence_path)
                @test read(target.uniprot_sequence_path, String) == ">P20963\nACDE\n"
                @test read(target.ensembl_protein_sequence_path, String) ==
                      ">ENSP000001\nACDE\n"
                @test calls[] == 1
            end
        end
    end

    @testset "unfetched override completes and reuses the real target stage" begin
        mktempdir() do root
            calls = Ref(0)
            resolver = override_resolver(calls)
            options = (; stage_options..., ensembl_protein_id = "ENSP000002", resolver)
            target = Iduna._resolve_target_stage("P20963", root; options...)
            @test target.ensembl_protein_id == "ENSP000002"
            @test target.ensembl_protein_sequence_path === nothing
            @test target.sequence_validated === nothing
            @test target.mapping_confirmed === nothing
            @test any(contains("override"), target.warnings)
            @test isfile(target.uniprot_sequence_path)
            @test isfile(joinpath(root, "sequences", "ensembl_proteins", "ENSP000001.fasta"))
            @test !ispath(joinpath(root, "sequences", "ensembl_proteins", "ENSP000002.fasta"))

            state_path = Iduna.Utils._stage_state_path(Iduna._target_stage_dir(root))
            state = JSON.parse(read(state_path, String))
            saved = JSON.parse(read(joinpath(root, "target.json"), String))
            @test state["status"] == "done"
            @test saved["ensembl_protein_sequence_path"] === nothing
            @test state["required_outputs"]["ensembl_sequence"] === nothing
            @test Set(keys(state["output_digests"])) ==
                  Set(["target.json", joinpath("sequences", "uniprot", "P20963.fasta")])

            cached = Iduna._resolve_target_stage("P20963", root; options...)
            @test calls[] == 1
            @test cached.ensembl_protein_id == target.ensembl_protein_id
            @test cached.ensembl_protein_sequence_path === nothing
            @test cached.warnings == target.warnings
            @test JSON.parse(read(state_path, String))["action"] == "reuse"

            write(target.uniprot_sequence_path, ">P20963\nACDF\n")
            rebuilt = Iduna._resolve_target_stage("P20963", root; options...)
            @test calls[] == 2
            @test read(rebuilt.uniprot_sequence_path, String) == ">P20963\nACDE\n"
            @test rebuilt.ensembl_protein_sequence_path === nothing
        end
    end

    @testset "fetched override remains an integrity-checked artifact" begin
        mktempdir() do root
            calls = Ref(0)
            resolver = override_resolver(calls; fetch_override = true)
            options = (; stage_options..., ensembl_protein_id = "ENSP000002", resolver)
            target = Iduna._resolve_target_stage("P20963", root; options...)
            protein_path = joinpath(root, "sequences", "ensembl_proteins", "ENSP000002.fasta")
            @test target.ensembl_protein_sequence_path == protein_path
            @test target.sequence_validated === nothing
            @test target.mapping_confirmed === nothing
            state_path = Iduna.Utils._stage_state_path(Iduna._target_stage_dir(root))
            state = JSON.parse(read(state_path, String))
            relative_protein = relpath(protein_path, root)
            @test state["required_outputs"]["ensembl_sequence"] == relative_protein
            @test state["output_digests"][relative_protein] ==
                  Iduna.Utils._file_sha256(protein_path)
            Iduna._resolve_target_stage("P20963", root; options...)
            @test calls[] == 1
            write(protein_path, ">ENSP000002\nACDF\n")
            rebuilt = Iduna._resolve_target_stage("P20963", root; options...)
            @test calls[] == 2
            @test read(rebuilt.ensembl_protein_sequence_path, String) ==
                  ">ENSP000002\nACDE\n"
            @test JSON.parse(read(state_path, String))["action"] == "rebuild"
        end
    end

    @testset "claimed sequence artifacts still must be complete" begin
        for (protein, fetch_override, damage) in (
            ("ENSP000001", false,
            dir -> rm(joinpath(dir, "ensembl_proteins", "ENSP000001.fasta"))),
            ("ENSP000002", true,
            dir -> write(joinpath(dir, "ensembl_proteins", "ENSP000002.fasta"), "invalid")),
            ("ENSP000002", false, dir -> rm(joinpath(dir, "uniprot", "P20963.fasta"))))
            mktempdir() do root
                calls = Ref(0)
                resolver = override_resolver(calls; fetch_override, damage)
                @test_throws "Resolved target sequence artifacts are missing." Iduna._resolve_target_stage(
                    "P20963", root; stage_options..., ensembl_protein_id = protein, resolver)
                @test calls[] == 1
                state_path = Iduna.Utils._stage_state_path(Iduna._target_stage_dir(root))
                @test JSON.parse(read(state_path, String))["status"] == "failed"
                @test !isfile(joinpath(root, "target.json"))
            end
        end
    end
end
