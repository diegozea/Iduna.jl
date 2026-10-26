let
    # Keep external-tool logs if a regression fails, including inside the full suite.
    root = mktempdir(; cleanup = false)
    @info "Local cluster-member regression artifacts." directory=root
    result = @testset "real cluster members and seed columns" begin
        expansion = Iduna.MSAExpansion
        msa = Iduna.ResultsValidation
        sequence = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALY" *
                   "LVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        # Keep changes internal: hmmalign --trim can remove unaligned hit termini.
        member1 = sequence[1:70] * "A" * sequence[72:end]
        member2 = sequence[1:40] * "WW" * sequence[41:end]
        representative_fasta = joinpath(root, "representative.fasta")
        members_fasta = joinpath(root, "members.fasta")
        write(representative_fasta, ">representative\n$(sequence)\n")
        write(members_fasta,
            ">representative\n$(sequence)\n>member1\n$(member1)\n>member2\n$(member2)\n")

        db = joinpath(root, "cluster_db")
        mmseqs = expansion.MMseqs2_jll.mmseqs()
        logs = joinpath(root, "database_logs")
        expansion._run_labeled(`$mmseqs createdb $representative_fasta $db`,
            "representative", logs)
        expansion._run_labeled(`$mmseqs createdb $members_fasta $(db * "_seq")`,
            "members", logs)
        expansion._run_labeled(
            `$mmseqs search $db $(db * "_seq") $(db * "_aln") $(joinpath(root, "search_tmp")) -a --threads 1`,
            "cluster_members", logs)

        # Two rows lack query residues 61–62: --symfrac 0.0 must retain those states.
        # Follow query residues; the separate physical all-gap columns may disappear.
        sparse_seed = sequence[1:60] * "--" * sequence[63:end]
        add_gap_columns(s) = s[1:40] * "--" * s[41:end]
        seed_codes = repeat("a", 40) * "bb" * repeat("c", 70)
        seed_path = joinpath(root, "seed.sto")
        write(seed_path,
            "# STOCKHOLM 1.0\n" *
            "#=GF SExonCodeMap \"a\"=>\"1_0\",\"b\"=>\"2_0\",\"c\"=>\"3_0\"\n" *
            "ENSG1 $(add_gap_columns(sequence))\n" *
            "seed_gap1 $(add_gap_columns(sparse_seed))\n" *
            "seed_gap2 $(add_gap_columns(sparse_seed))\n" *
            "#=GC SExonCode $(seed_codes)\n//\n")
        target = Iduna.ResolvedTarget(;
            input_id = "ENST1", input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG1", transcript_id = "ENST1")
        seed = Iduna.SeedSelection(; pid = 10.0, epli = 100.0,
            stockholm_path = seed_path, summary_path = joinpath(root, "summary.csv"))
        expanded = expansion.expand_msa(target, seed, root;
            mmseqs_db = db, mmseqs_threads = 1, centroids = true)

        full = msa.read_file(expanded.full_stockholm, msa.Stockholm;
            keepinserts = true, deletefullgaps = false)
        match_only = msa.read_file(expanded.match_stockholm, msa.Stockholm;
            keepinserts = true, deletefullgaps = false)
        expected_sequences = Dict("ENSG1" => sequence,
            "seed_gap1" => replace(sparse_seed, "-" => ""),
            "seed_gap2" => replace(sparse_seed, "-" => ""),
            "representative" => sequence, "member1" => member1, "member2" => member2)
        @test expanded.status === :ok
        @test expanded.n_hits == 3
        @test expanded.n_new_hits == 3
        names = Set(msa.sequencenames(full))
        @test names == Set(keys(expected_sequences))
        @test Set(msa.sequencenames(match_only)) == Set(keys(expected_sequences))
        for (name, expected) in expected_sequences
            @test name in names &&
                  replace(msa.stringsequence(full, name), "-" => "") == expected
        end

        # Labels follow query residues; added hit insertions remain unassigned.
        query = msa.stringsequence(full, "ENSG1")
        labels = collect(Iduna.Utils.s_exon_codes(full))
        query_columns = findall(!=('-'), collect(query))
        @test length(labels) == msa.ncolumns(full)
        @test String(labels[query_columns]) == repeat("a", 40) * repeat("c", 70)
        @test all(==('.'), labels[findall(==('-'), collect(query))])
        @test msa.ncolumns(match_only) == 110
        @test msa.stringsequence(match_only, "ENSG1") == sequence
        @test Iduna.Utils.s_exon_codes(match_only) == repeat("a", 40) * repeat("c", 70)
        @test Iduna.Utils.s_exon_code_map(full) ==
              Dict('a' => "1_0", 'b' => "2_0", 'c' => "3_0")

        centroids = msa.read_file(
            joinpath(expanded.run_dir, "centroid_msa",
                "ENST1_centroids_full.sto"),
            msa.Stockholm; keepinserts = true)
        @test Set(msa.sequencenames(centroids)) ==
              Set(["ENSG1", "seed_gap1", "seed_gap2", "representative"])
    end
    counts = Test.get_test_counts(result)
    if counts.fails == 0 && counts.errors == 0
        rm(root; recursive = true)
    end
end
