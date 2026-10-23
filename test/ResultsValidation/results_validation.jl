@testset "ResultsValidation" begin
    mktempdir() do tmp
        sto = joinpath(tmp, "tiny.sto")
        write(sto, """
        # STOCKHOLM 1.0
        seq1 ACDE
        seq2 ACDF
        //
        """)
        msa = Iduna.ResultsValidation.load_msa(sto)
        @test Iduna.ResultsValidation.nsequences(msa) == 2
        @test_throws ErrorException Iduna.ResultsValidation._resolve_query_name(
            msa, "missing_gene", "missing_tx")
        stats = Iduna.ResultsValidation.alignment_stats(sto)
        @test stats.n_sequences == 2
        @test stats.n_columns == 4

        a3m = joinpath(tmp, "tiny.a3m")
        write(a3m, ">seq1\nACde\n>seq2\nACfg\n")
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_msa(a3m)) == 2

        fasta = joinpath(tmp, "tiny.fasta")
        write(fasta, ">seq1\nACDE\n>seq2\nACDF\n")
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_msa(fasta)) == 2
        @test_throws ErrorException Iduna.ResultsValidation.load_msa(
            joinpath(tmp, "unknown.txt"))

        seed = Iduna.SeedSelection(;
            pid = 10.0,
            epli = 100.0,
            stockholm_path = sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        expansion = Iduna.ExpansionResult(;
            run_dir = tmp,
            seed_stockholm = sto,
            hits_fasta = fasta,
            full_stockholm = sto,
            match_stockholm = sto,
            a3m_path = a3m,
            db_dir = tmp,
            hmm_dir = tmp,
            logs_dir = tmp
        )
        uniprot = joinpath(tmp, "uniprot.fasta")
        write(uniprot, ">P20963\nACDF\n")
        target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001",
            uniprot_sequence_path = uniprot
        )
        validation = Iduna.ResultsValidation.validate_results(target, seed, expansion, tmp)
        @test validation.status === :warn
        @test validation.query_name == "seq1"
        @test validation.aln_mismatches == 1
        @test validation.aln_insertions == 0
        @test validation.aln_deletions == 0
        @test occursin("substitutions", only(validation.warnings))
        @test isfile(validation.stats_path)
        @test isfile(validation.query_vs_uniprot_path)
        stats_row = only(collect(Iduna.ResultsValidation.CSV.File(validation.stats_path)))
        @test stats_row.seed_path == "tiny.sto"
        @test stats_row.expanded_path == "tiny.sto"
        @test stats_row.query_vs_uniprot_path ==
              joinpath("validation", "pid_10.00", "query_vs_uniprot_alignment.txt")

        original_stats = read(validation.stats_path, String)
        write(validation.stats_path, replace(original_stats, "seq1" => "seqX"))
        rebuilt = Iduna.ResultsValidation.validate_results(target, seed, expansion, tmp)
        @test rebuilt.query_name == "seq1"
        @test read(rebuilt.stats_path, String) == original_stats
        @test Iduna.Utils._read_stage_state(dirname(rebuilt.stats_path))["action"] ==
              "rebuild"
        write(rebuilt.query_vs_uniprot_path, "")
        rebuilt = Iduna.ResultsValidation.validate_results(target, seed, expansion, tmp)
        @test filesize(rebuilt.query_vs_uniprot_path) > 0
        @test Iduna.Utils._read_stage_state(dirname(rebuilt.stats_path))["action"] ==
              "rebuild"
        unchanged = Iduna.ResultsValidation.validate_results(target, seed, expansion, tmp)
        @test unchanged.query_name == "seq1"
        @test Iduna.Utils._read_stage_state(dirname(unchanged.stats_path))["action"] ==
              "reuse"

        paths = Iduna.ResultsValidation._validation_input_paths(target, seed, expansion, tmp)
        outputs = Iduna.ResultsValidation._validation_outputs(tmp, seed, paths)
        identity = Iduna.ResultsValidation._validation_identity(target, seed, expansion, paths)
        write(outputs.query_vs_uniprot_path, "")
        @test_throws ErrorException Iduna.ResultsValidation._write_validation_state(
            tmp, seed, :done, identity, outputs; action = :run)
        write(outputs.stats_path, "query_name\n")
        @test_throws ErrorException Iduna.ResultsValidation._write_validation_state(
            tmp, seed, :done, identity, outputs; action = :run)

        indel_uniprot = joinpath(tmp, "uniprot_indel.fasta")
        write(indel_uniprot, ">P20963\nACD\n")
        indel_target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001",
            uniprot_sequence_path = indel_uniprot
        )
        indel_validation = Iduna.ResultsValidation.validate_results(
            indel_target, seed, expansion, joinpath(tmp, "indel"))
        @test indel_validation.status === :warn
        @test occursin("indels relative to the UniProt sequence",
            only(indel_validation.warnings))
        indel_row = (;
            aln_identical = false,
            aln_insertions = 1,
            aln_deletions = 0)
        @test occursin("Seed query has indels",
            only(Iduna.ResultsValidation._cached_alignment_warnings(indel_row, nothing)))
        identical_row = (;
            aln_identical = true,
            aln_insertions = 0,
            aln_deletions = 0)
        @test isempty(
            Iduna.ResultsValidation._cached_alignment_warnings(identical_row, nothing))
        @test isempty(Iduna.ResultsValidation._alignment_warnings(nothing,
            (; insertions = 0, deletions = 0, identical = true)))

        relative_seed = Iduna.SeedSelection(;
            pid = 10.0,
            epli = 100.0,
            stockholm_path = "tiny.sto",
            summary_path = "seed_summary.csv",
            workdir = tmp)
        relative_expansion = Iduna.ExpansionResult(;
            run_dir = ".",
            seed_stockholm = "tiny.sto",
            hits_fasta = "tiny.fasta",
            full_stockholm = "tiny.sto",
            match_stockholm = "tiny.sto",
            a3m_path = "tiny.a3m",
            db_dir = ".",
            hmm_dir = ".",
            logs_dir = ".",
            workdir = tmp)
        relative_target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001",
            uniprot_sequence_path = "uniprot.fasta",
            workdir = tmp)
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_seed_msa(relative_seed)) == 2
        @test Iduna.ResultsValidation.nsequences(
            Iduna.ResultsValidation.load_expanded_msa(relative_expansion)) == 2
        relative_validation = Iduna.ResultsValidation.validate_results(
            relative_target, relative_seed, relative_expansion, tmp)
        @test relative_validation.status === :warn
        @test relative_validation.query_name == "seq1"

        seed_only_validation = Iduna.ResultsValidation.validate_results(
            target, seed, nothing, joinpath(tmp, "seed_only"))
        @test seed_only_validation.status === :warn
        @test seed_only_validation.query_name == "seq1"
        @test seed_only_validation.seed_nseq == 2
        @test seed_only_validation.expanded_nseq === nothing
        @test seed_only_validation.expanded_ncol === nothing
        @test seed_only_validation.expanded_clusters62 === nothing
        @test seed_only_validation.expanded_neff80 === nothing
        @test occursin("Seed query has substitutions", only(seed_only_validation.warnings))
        @test isfile(seed_only_validation.stats_path)
        @test isfile(seed_only_validation.query_vs_uniprot_path)

        transcript_target = Iduna.ResolvedTarget(;
            input_id = "ENST000001",
            input_kind = :ensembl_transcript,
            uniprot_id = "P20963",
            ensembl_gene_id = "seq1",
            transcript_id = "ENST000001"
        )
        transcript_validation = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, expansion, joinpath(tmp, "transcript"))
        @test transcript_validation.status === :ok
        @test isfile(transcript_validation.stats_path)
        overwritten_transcript_validation = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, expansion, joinpath(tmp, "transcript");
            overwrite = true)
        @test overwritten_transcript_validation.status === :ok
        cached_transcript_validation = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, expansion, joinpath(tmp, "transcript"))
        @test cached_transcript_validation.status === :ok
        skipped_expansion = Iduna.ExpansionResult(;
            run_dir = expansion.run_dir,
            seed_stockholm = expansion.seed_stockholm,
            seed_fasta = expansion.seed_fasta,
            hits_fasta = expansion.hits_fasta,
            full_stockholm = expansion.full_stockholm,
            match_stockholm = expansion.match_stockholm,
            a3m_path = expansion.a3m_path,
            s_exon_blocks_tsv = expansion.s_exon_blocks_tsv,
            db_dir = expansion.db_dir,
            hmm_dir = expansion.hmm_dir,
            logs_dir = expansion.logs_dir,
            n_hits = expansion.n_hits,
            n_new_hits = expansion.n_new_hits,
            status = :skipped,
            workdir = expansion.workdir)
        cached_after_expansion_reuse = Iduna.ResultsValidation.validate_results(
            transcript_target, seed, skipped_expansion, joinpath(tmp, "transcript"))
        @test cached_after_expansion_reuse.status === :ok
        validation_state = Iduna.Utils._read_stage_state(
            Iduna.ResultsValidation._validation_dir(joinpath(tmp, "transcript"), seed))
        @test validation_state["action"] == "reuse"

        missing_seed = Iduna.SeedSelection(;
            pid = 20.0,
            epli = 100.0,
            stockholm_path = joinpath(tmp, "missing_seed.sto"),
            summary_path = joinpath(tmp, "missing_summary.csv"))
        failure_workdir = joinpath(tmp, "validation_failure")
        @test_throws ErrorException Iduna.ResultsValidation.validate_results(
            transcript_target, missing_seed, nothing, failure_workdir)
        failed_state = Iduna.Utils._read_stage_state(
            Iduna.ResultsValidation._validation_dir(failure_workdir, missing_seed))
        @test failed_state["status"] == "failed"
        @test failed_state["exception"]["type"] == "ErrorException"
    end
end

@testset "validation certifies only the requested query" begin
    mktempdir() do tmp
        uniprot = joinpath(tmp, "uniprot.fasta")
        write(uniprot, ">reference\nACDE\n")
        for (label, rows, succeeds) in (
            ("absent", "wrong ACDE\nother ACDF\n", false),
            ("gene", "ENSG1 ACDE\nother ACDF\n", true),
            ("transcript", "ENST1 ACDE\nother ACDF\n", true),
            ("variant", "1|ENSG1.2 ACDE\nother ACDF\n", true),
            ("ambiguous", "1|ENSG1.2 ACDE\n2|ENSG1.3 ACDE\n", false),
            ("two_ids", "ENSG1 ACDE\nENST1 ACDE\n", false))
            sto = joinpath(tmp, "$(label).sto")
            write(sto, "# STOCKHOLM 1.0\n$(rows)//\n")
            seed = Iduna.SeedSelection(pid = 10.0, epli = 100.0,
                stockholm_path = sto, summary_path = joinpath(tmp, "summary.csv"))
            for reference in (nothing, uniprot)
                target = Iduna.ResolvedTarget(input_id = "ENST1",
                    input_kind = :ensembl_transcript, ensembl_gene_id = "ENSG1",
                    transcript_id = "ENST1", uniprot_sequence_path = reference)
                workdir = joinpath(tmp, label, reference === nothing ? "no_ref" : "ref")
                if succeeds
                    result = Iduna.ResultsValidation.validate_results(target, seed, nothing, workdir)
                    @test result.status === :ok
                    @test result.aln_identical === (reference === nothing ? nothing : true)
                else
                    @test_throws ErrorException Iduna.ResultsValidation.validate_results(
                        target, seed, nothing, workdir)
                    state = Iduna.Utils._read_stage_state(
                        Iduna.ResultsValidation._validation_dir(workdir, seed))
                    @test state["status"] == "failed"
                    @test occursin("query", state["exception"]["message"])
                end
            end
        end
        seed_path = joinpath(tmp, "gene.sto")
        seed = Iduna.SeedSelection(pid = 10.0, epli = 100.0,
            stockholm_path = seed_path, summary_path = joinpath(tmp, "summary.csv"))
        expansion = Iduna.ExpansionResult(run_dir = tmp, seed_stockholm = seed_path,
            hits_fasta = uniprot, full_stockholm = joinpath(tmp, "absent.sto"),
            match_stockholm = joinpath(tmp, "absent.sto"), a3m_path = uniprot,
            db_dir = tmp, hmm_dir = tmp, logs_dir = tmp)
        target = Iduna.ResolvedTarget(input_id = "ENST1", input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG1", transcript_id = "ENST1", uniprot_sequence_path = uniprot)
        @test_throws ErrorException Iduna.ResultsValidation.validate_results(
            target, seed, expansion, joinpath(tmp, "bad_expansion"))

        # Saved directional columns must use the query-relative convention too.
        write(seed_path, "# STOCKHOLM 1.0\nENSG1 ACDEF\nother ACDEF\n//\n")
        result = Iduna.ResultsValidation.validate_results(target, seed, nothing, joinpath(tmp, "indel"))
        @test (result.aln_insertions, result.aln_deletions) == (1, 0)
        row = only(Iduna.ResultsValidation.CSV.File(result.stats_path))
        @test (row.aln_insertions, row.aln_deletions) == (1, 0)
        @test occursin("insertions=1", read(result.query_vs_uniprot_path, String))

        # A pre-fix completion record must not reuse reversed saved diagnostics.
        workdir = joinpath(tmp, "indel")
        stage_dir = Iduna.ResultsValidation._validation_dir(workdir, seed)
        state = Iduna.Utils._read_stage_state(stage_dir)
        delete!(state["identity"], "validation_contract")
        state["identity_hash"] = Iduna.Utils._identity_hash(state["identity"])
        Iduna.Utils.write_json(Iduna.Utils._stage_state_path(stage_dir), state)
        rebuilt = Iduna.ResultsValidation.validate_results(target, seed, nothing, workdir)
        @test (rebuilt.aln_insertions, rebuilt.aln_deletions) == (1, 0)
        @test Iduna.Utils._read_stage_state(stage_dir)["action"] == "rebuild"
    end
end
