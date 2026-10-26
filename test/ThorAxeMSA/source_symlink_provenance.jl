function test_thoraxe_source_symlinks(root, package, python)
    msa = Iduna.ThorAxeMSA
    utils = Iduna.Utils
    @testset "ThorAxe linked scientific source and public cache" begin
        score() = parse(Float64,
            strip(read(
                `$python -B -c "from thoraxe.subexons.alignment import _percent_identity_without_gaps; print(_percent_identity_without_gaps('AAAA', 'AAAA'))"`,
                String)))
        source = joinpath(root, "cached_input")
        ensembl = mkpath(joinpath(source, "Ensembl"))
        for file in msa._REQUIRED_ENSEMBL_FILES
            write(joinpath(ensembl, file), "$(file)\n")
        end
        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1", input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1", transcript_id = "ENST00000000001.1")
        calls = Ref(0)
        function scientific_fixture(input, run_root; kwargs...)
            calls[] += 1
            # Only replace the biological CLI; use its real imported function to
            # make changes in the copied scientific source observable in the MSA.
            sequence = score() == 100.0 ? "AA" : "AC"
            thoraxe = joinpath(run_root, "thoraxe")
            mkpath(joinpath(thoraxe, "msa"))
            mkpath(joinpath(thoraxe, "phylosofs"))
            write(joinpath(thoraxe, "path_table.csv"),
                "TranscriptIDCluster,Path\nENST00000000001,start/1_0/stop\n")
            write(joinpath(thoraxe, "s_exon_table.csv"),
                "GeneID,Species,TranscriptIDCluster,S_exonID,S_exon_Sequence\n" *
                "ENSG00000000001,homo_sapiens,ENST00000000001,1_0,$sequence\n" *
                "ORTHO1,mus_musculus,ORTHO1,1_0,AB\n")
            write(joinpath(thoraxe, "msa", "msa_s_exon_1_0.fasta"),
                ">ENSG00000000001\n$sequence\n>ORTHO1\nAB\n")
            write(joinpath(thoraxe, "phylosofs", "s_exons.tsv"), "1_0\ta\n")
            write(joinpath(thoraxe, "phylosofs", "transcripts.pir"),
                ">P1;ENSG00000000001 ENST00000000001 a\naa\n$sequence*\n")
            return thoraxe
        end
        runner = (args...; kwargs...) -> msa._run_thoraxe_msa_stage!(
            args...; kwargs..., thoraxe_fn = scientific_fixture)
        workdir = joinpath(root, "linked_source_work")
        function build(path = workdir; overwrite = false)
            result = msa.build_thoraxe_msa(target, path;
                pid_thresholds = [10.0], specieslist = "all",
                cached_thoraxe_input_dir = source, specieslist_filter = false,
                biomart_datasets_filter = false, pid_sample_count = 0, pid_sample_seed = 7,
                overwrite, _stage_runner = runner)
            @test result.status === :ok
            return read(only(result.baseline_fastas), String)
        end
        state() = utils._read_stage_state(msa._thoraxe_msa_stage_dir(workdir))

        @test score() == 100.0
        plain = msa._thoraxe_stage_provenance(0)
        @test msa._known_thoraxe_provenance(plain, 0)
        before = build()
        @test build() == before
        @test calls[] == 1
        @test state()["action"] == "reuse"

        # A symlink to a readable source file is supported and hashes its bytes.
        module_path = joinpath(package, "subexons", "alignment.py")
        file_source = joinpath(root, "alignment.py")
        mv(module_path, file_source)
        original_source = read(file_source, String)
        old_return = "return 100.0 * (identical / aln_len)"
        @test count(old_return, original_source) == 1
        changed_source = replace(original_source, old_return => "return 75.0 * (identical / aln_len)")
        symlink(file_source, module_path)
        @test msa._thoraxe_stage_provenance(0) == plain
        @test build() == before
        @test calls[] == 1
        write(file_source, changed_source)
        @test score() == 75.0
        changed = msa._thoraxe_stage_provenance(0)
        @test msa._known_thoraxe_provenance(changed, 0)
        @test changed.python != plain.python
        after = build()
        @test after != before
        @test calls[] == 2
        @test build() == after
        @test calls[] == 2

        # Retarget a directory link between immutable source snapshots. Each
        # scientific observation uses a fresh Python process, with no hot reload.
        rm(module_path)
        write(module_path, original_source)
        subpackage = joinpath(package, "subexons")
        first_snapshot = joinpath(root, "subexons_first")
        second_snapshot = joinpath(root, "subexons_second")
        mv(subpackage, first_snapshot)
        cp(first_snapshot, second_snapshot)
        write(joinpath(second_snapshot, "alignment.py"), changed_source)
        snapshots = [utils._file_sha256(joinpath(path, "alignment.py"))
                     for path in (first_snapshot, second_snapshot)]
        symlink(first_snapshot, subpackage)
        for (snapshot, expected_score, expected_msa) in ((first_snapshot, 100.0, before), (
            second_snapshot, 75.0, after))
            rm(subpackage)
            symlink(snapshot, subpackage)
            @test score() == expected_score
            provenance = msa._thoraxe_stage_provenance(0)
            @test provenance.python === nothing
            @test !msa._known_thoraxe_provenance(provenance, 0)
            previous_calls = calls[]
            @test build() == expected_msa
            @test state()["action"] == "rebuild"
            @test state()["status"] == "done"
            @test build() == expected_msa
            @test calls[] == previous_calls + 2
            @test build(joinpath(root, "fresh_$(Int(expected_score))")) == expected_msa
        end
        @test snapshots == [utils._file_sha256(joinpath(path, "alignment.py"))
               for path in (first_snapshot, second_snapshot)]

        rm(subpackage)
        mv(first_snapshot, subpackage)
        @test msa._known_thoraxe_provenance(msa._thoraxe_stage_provenance(0), 0)
        # Unreadable link targets must not disappear from the inventory either.
        for destination in (joinpath(root, "missing_source"), joinpath(package, "cycle"))
            link = joinpath(package, "cycle")
            symlink(destination, link)
            @test msa._transcript_query_tool_identity() === nothing
            rm(link)
        end
    end
end
