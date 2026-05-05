@testset "MSAExpansion" begin
    mktempdir() do tmp
        source = joinpath(tmp, "in.sto")
        dest = joinpath(tmp, "out.sto")
        write(source, "seq1 AC\n")
        Iduna.MSAExpansion.prepare_stockholm_for_mmseqs(source, dest)
        text = read(dest, String)
        @test startswith(text, "# STOCKHOLM 1.0")
        @test endswith(chomp(text), "//")

        fragmented = joinpath(tmp, "fragmented.sto")
        write(fragmented, """
        # STOCKHOLM 1.0
        # source comment
        #=GF ID family1
        #=GF DE first description
        #=GS seq1 AC accession1
        seq1 AC
        #=GR seq1 PP 99
        #=GC RF xx
        #=GF
        #=GS seq1
        seq1 DE
        #=GR seq1 PP **
        #=GC RF yy
        //
        """)
        Iduna.MSAExpansion.normalize_stockholm_annotations!(fragmented)
        normalized = read(fragmented, String)
        @test occursin("# source comment", normalized)
        @test occursin("#=GF ID family1", normalized)
        @test occursin("#=GF DE first description", normalized)
        @test occursin("#=GS seq1 AC accession1", normalized)
        @test occursin("seq1\tACDE", normalized)
        @test occursin("#=GC RF xxyy", normalized)
        @test occursin("#=GR seq1 PP 99**", normalized)

        hits_tsv = joinpath(tmp, "hits.tsv")
        write(hits_tsv, "query\tseed one\tACD-\nquery\thit one\tACDF\n")
        all_hits, filtered_hits = Iduna.MSAExpansion._collect_hits(hits_tsv, Set(["seed"]))
        @test all_hits == [("seed", "ACD"), ("hit", "ACDF")]
        @test filtered_hits == [("hit", "ACDF")]

        empty_hits = joinpath(tmp, "empty_hits.fasta")
        touch(empty_hits)
        @test Iduna.MSAExpansion._cached_hit_counts(empty_hits, Set(["seed"])) ==
              (n_hits = 0, n_new_hits = 0)

        reorder_sto = joinpath(tmp, "reorder.sto")
        write(reorder_sto, "# STOCKHOLM 1.0\nb AC\na AC\nc AC\n//\n")
        reordered = Iduna.MSAExpansion._reorder_alignment(
            Iduna.MSAExpansion.read_file(reorder_sto, Iduna.MSAExpansion.Stockholm),
            ["a", "missing"])
        @test String.(Iduna.MSAExpansion.sequencenames(reordered)) == ["a", "b", "c"]

        gene_id = "ENSG00000198821"
        transcript_id = "ENST00000362089.10"
        seed_sto = joinpath(tmp, "seed.sto")
        write(seed_sto, "# STOCKHOLM 1.0\nseed ACDE\n//\n")
        seed = Iduna.SeedSelection(;
            pid = 10.0,
            median_identity = 100.0,
            mean_identity = 100.0,
            stockholm_path = seed_sto,
            summary_path = joinpath(tmp, "seed_summary.csv")
        )
        target = Iduna.ResolvedTarget(;
            input_id = transcript_id,
            input_kind = :ensembl_transcript,
            ensembl_gene_id = gene_id,
            transcript_id
        )

        db = joinpath(tmp, "mock_mmseqs_db")
        touch("$(db).dbtype")
        touch("$(db)_aln.dbtype")
        touch("$(db)_seq.dbtype")

        cached_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "expanded_msa")
        mkpath(cached_dir)
        write(joinpath(cached_dir, "$(transcript_id)_full.sto"), "cached full\n")
        write(joinpath(cached_dir, "$(transcript_id)_matchonly.sto"), "cached match\n")
        write(joinpath(cached_dir, "$(transcript_id)_expanded.a3m"), "cached a3m\n")
        hits_fasta = joinpath(cached_dir, "$(transcript_id)_hits_raw.fasta")
        write(hits_fasta, ">seed one\nACDE\n>hit one\nACDF\n>hit_two\nACDG\n")

        cached = Iduna.MSAExpansion.expand_msa(target, seed, tmp; mmseqs_db = db)
        hits_msa = Iduna.MSAExpansion.read_file(hits_fasta, Iduna.MSAExpansion.FASTA)
        @test cached.status === :skipped
        @test cached.n_hits == Iduna.MSAExpansion.nsequences(hits_msa)
        @test cached.n_hits == 3
        @test cached.n_new_hits == 2
    end
end
