mktempdir() do tmp
    thoraxe = joinpath(tmp, "thoraxe")
    msa_dir = joinpath(thoraxe, "msa")
    mkpath(msa_dir)
    write(joinpath(thoraxe, "path_table.csv"),
        "TranscriptIDCluster,Path\nENST00000000001/ENST00000000002,start/1_0/2_0/stop\n")
    write(joinpath(thoraxe, "s_exon_table.csv"),
        "GeneID,Species\nENSG00000000001,homo_sapiens\nORTHO1,mus_musculus\n")
    write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
        ">ENSG00000000001\nA-A.\n>ORTHO1\nABCD\n")
    write(joinpath(msa_dir, "msa_s_exon_2_0.fasta"),
        ">ENSG00000000001\nCC\n>ORTHO1\nCD\n")
    phylosofs_dir = joinpath(thoraxe, "phylosofs")
    mkpath(phylosofs_dir)
    write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\ta\n2_0\t=\n")
    write(joinpath(phylosofs_dir, "transcripts.pir"), """
    >P1;ENSG00000000001 ENST00000000001/ENST00000000002 a=
    aa==
    AACC*
    """)

    msa,
    species = Iduna.ThorAxeMSA.assemble_transcript_msa(
        thoraxe, "ENSG00000000001.1", "ENST00000000002.1")
    @test Iduna.ThorAxeMSA.nsequences(msa) == 2
    name = first(Iduna.ThorAxeMSA.sequencenames(msa))
    @test length(Iduna.ThorAxeMSA.stringsequence(msa, name)) == 6
    @test Iduna.Utils.s_exon_codes(msa) == "aaaa=="
    @test Iduna.ThorAxeMSA._s_exon_symbol_for_reference_residue('-', 'a') == 'a'
    @test Iduna.ThorAxeMSA._s_exon_symbol_for_reference_residue('.', 'a') == '.'
    consumed_symbol,
    next_symbol_state = Iduna.ThorAxeMSA._consume_phylosofs_symbol(
        iterate("a"), "a", "1_0", Dict('a' => "1_0"), "ENST00000000001")
    @test consumed_symbol == 'a'
    @test next_symbol_state === nothing
    @test_throws ErrorException Iduna.ThorAxeMSA._consume_phylosofs_symbol(
        nothing, "a", "1_0", Dict('a' => "1_0"), "ENST00000000001")
    @test_throws ErrorException Iduna.ThorAxeMSA._consume_phylosofs_symbol(
        iterate("a"), "a", "2_0", Dict('a' => "1_0"), "ENST00000000001")
    @test Iduna.Utils.s_exon_code_map(msa) == Dict('a' => "1_0", '=' => "2_0")
    @test species == ["homo_sapiens", "mus_musculus"]

    uniprot = joinpath(tmp, "uniprot.fasta")
    target = Iduna.ResolvedTarget(;
        input_id = "P20963",
        input_kind = :uniprot,
        uniprot_id = "P20963",
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1",
        uniprot_sequence_path = uniprot)
    write(uniprot, ">P20963\nAACC\n")
    @test isempty(Iduna.ThorAxeMSA._validate_transcript_translation(target, msa))
    relative_target = Iduna.ResolvedTarget(;
        input_id = "P20963",
        input_kind = :uniprot,
        uniprot_id = "P20963",
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1",
        uniprot_sequence_path = "uniprot.fasta",
        workdir = tmp)
    @test isempty(Iduna.ThorAxeMSA._validate_transcript_translation(
        relative_target, msa; workdir = joinpath(tmp, "new_run")))

    write(uniprot, ">P20963\nAAAC\n")
    warnings = Iduna.ThorAxeMSA._validate_transcript_translation(target, msa)
    @test length(warnings) == 1
    @test occursin("substitutions", only(warnings))
    candidate_warning = Iduna.ThorAxeMSA._candidate_msa0_validation(target, msa, 80.0, tmp)
    @test candidate_warning.eligible === true
    @test candidate_warning.status == "warning"
    @test occursin("substitutions", candidate_warning.issue)

    write(uniprot, ">P20963\nAAC\n")
    @test_throws ErrorException Iduna.ThorAxeMSA._validate_transcript_translation(
        target, msa)
    candidate_invalid = Iduna.ThorAxeMSA._candidate_msa0_validation(target, msa, 80.0, tmp)
    @test candidate_invalid.eligible === false
    @test candidate_invalid.status == "invalid_msa0"
    @test occursin("indels", candidate_invalid.issue)

    missing_uniprot = Iduna.ResolvedTarget(;
        input_id = "P20963",
        input_kind = :uniprot,
        uniprot_id = "P20963",
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1",
        uniprot_sequence_path = joinpath(tmp, "missing.fasta"))
    missing_warnings = Iduna.ThorAxeMSA._validate_transcript_translation(
        missing_uniprot, msa)
    @test length(missing_warnings) == 1
    @test occursin("UniProt sequence file is missing", only(missing_warnings))

    no_uniprot = Iduna.ResolvedTarget(;
        input_id = "ENST00000000001.1",
        input_kind = :ensembl_transcript,
        ensembl_gene_id = "ENSG00000000001.1",
        transcript_id = "ENST00000000001.1")
    @test isempty(Iduna.ThorAxeMSA._validate_transcript_translation(no_uniprot, msa))

    missing_reference = Iduna.ResolvedTarget(;
        input_id = "P20963",
        input_kind = :uniprot,
        uniprot_id = "P20963",
        ensembl_gene_id = "ENSG_MISSING",
        transcript_id = "ENST_MISSING",
        uniprot_sequence_path = uniprot)
    @test_throws ErrorException Iduna.ThorAxeMSA._validate_transcript_translation(
        missing_reference, msa)
end

mktempdir() do tmp
    thoraxe = joinpath(tmp, "thoraxe")
    msa_dir = joinpath(thoraxe, "msa")
    mkpath(msa_dir)
    write(joinpath(thoraxe, "path_table.csv"),
        "TranscriptIDCluster,Path\nENST00000000001,start/0_0/1_0/stop\n")
    write(joinpath(thoraxe, "s_exon_table.csv"),
        "GeneID,Species,TranscriptIDCluster,S_exonID,S_exon_Sequence\n" *
        "ENSG00000000001,homo_sapiens,ENST00000000001,0_0,MM\n" *
        "ENSG00000000001,homo_sapiens,ENST00000000001,1_0,AA\n" *
        "ORTHO1,mus_musculus,ORTHO1,1_0,AB\n")
    write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
        ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")
    phylosofs_dir = joinpath(thoraxe, "phylosofs")
    mkpath(phylosofs_dir)
    write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\ta\n")
    write(joinpath(phylosofs_dir, "transcripts.pir"), """
    >P1;ENSG00000000001 ENST00000000001 za
    zzaa
    MMAA*
    """)

    fallback_msa = Iduna.ThorAxeMSA._transcript_exon_msa(
        thoraxe, "0_0", "ENSG00000000001.1", "ENST00000000001.1")
    @test Iduna.ThorAxeMSA.nsequences(fallback_msa) == 1
    fallback_name = only(Iduna.ThorAxeMSA.sequencenames(fallback_msa))
    @test Iduna.ThorAxeMSA.stringsequence(fallback_msa, fallback_name) == "MM"

    msa,
    species = Iduna.ThorAxeMSA.assemble_transcript_msa(
        thoraxe, "ENSG00000000001.1", "ENST00000000001.1")
    reference = Iduna.ThorAxeMSA.resolve_sequence_name(msa, "ENSG00000000001")
    @test Iduna.ThorAxeMSA.stringsequence(msa, reference) == "MMAA"
    codes = Iduna.Utils.s_exon_codes(msa)
    @test length(codes) == 4
    @test codes[3:4] == "aa"
    @test codes[1] == codes[2]
    @test codes[1] != 'a'
    @test Iduna.Utils.s_exon_code_map(msa)[codes[1]] == "0_0"
    @test Iduna.Utils.s_exon_code_map(msa)['a'] == "1_0"
    @test species == ["homo_sapiens", "mus_musculus"]
end

mktempdir() do tmp
    thoraxe = joinpath(tmp, "thoraxe")
    msa_dir = joinpath(thoraxe, "msa")
    mkpath(msa_dir)
    write(joinpath(thoraxe, "path_table.csv"),
        "TranscriptIDCluster,Path\nENST00000000001,start/1_0/0_1/stop\n")
    write(joinpath(thoraxe, "s_exon_table.csv"),
        "GeneID,Species,TranscriptIDCluster,S_exonID,S_exon_Sequence\n" *
        "ENSG00000000001,homo_sapiens,ENST00000000001,1_0,AA\n" *
        "ENSG00000000001,homo_sapiens,ENST00000000001,0_1,\n")
    write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
        ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")
    phylosofs_dir = joinpath(thoraxe, "phylosofs")
    mkpath(phylosofs_dir)
    write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\ta\n0_1\tb\n")
    write(joinpath(phylosofs_dir, "transcripts.pir"), """
    >P1;ENSG00000000001 ENST00000000001 ab
    aa
    AA*
    """)

    @test Iduna.ThorAxeMSA._s_exon_sequence_value(missing) == ""
    @test Iduna.ThorAxeMSA._s_exon_sequence_value("") == ""
    @test Iduna.ThorAxeMSA._s_exon_sequence_value("A*") == "A"
    msa,
    species = Iduna.ThorAxeMSA.assemble_transcript_msa(
        thoraxe, "ENSG00000000001.1", "ENST00000000001.1")
    reference = Iduna.ThorAxeMSA.resolve_sequence_name(msa, "ENSG00000000001")
    @test Iduna.ThorAxeMSA.stringsequence(msa, reference) == "AA"
    @test Iduna.Utils.s_exon_codes(msa) == "aa"
    @test Iduna.Utils.s_exon_code_map(msa)['a'] == "1_0"
    @test Iduna.Utils.s_exon_code_map(msa)['b'] == "0_1"
    @test species == ["homo_sapiens", "unknown"]
end

mktempdir() do tmp
    thoraxe = joinpath(tmp, "thoraxe")
    msa_dir = joinpath(thoraxe, "msa")
    mkpath(msa_dir)
    write(joinpath(thoraxe, "path_table.csv"),
        "TranscriptIDCluster,Path\nENST00000000001,start/1_0/stop\n")
    write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
        ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")

    msa,
    species = Iduna.ThorAxeMSA.assemble_transcript_msa(
        thoraxe, "ENSG00000000001.1", "ENST00000000001.1")
    @test Iduna.ThorAxeMSA.nsequences(msa) == 2
    @test !Iduna.Utils.has_s_exon_annotations(msa)
    @test species == ["unknown", "unknown"]
end

mktempdir() do tmp
    summary = joinpath(tmp, "candidate_summary.csv")
    write(summary,
        "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
        "30.0,1,false,90.0,20,pid30.sto,pid30.fa\n" *
        "10.0,2,true,70.0,10,pid10.sto,pid10.fa\n" *
        "80.0,3,true,60.0,30,pid80.sto,pid80.fa\n")
    seed = Iduna.ThorAxeMSA.select_best_seed(summary)
    @test seed.pid == 10.0
    @test seed.epli == 70.0
    @test seed.workdir == abspath(tmp)

    write(summary,
        "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
        "10.0,1,true,70.0,10,pid10.sto,pid10.fa\n" *
        "80.0,2,true,70.0,20,pid80.sto,pid80.fa\n")
    larger_msa_seed = Iduna.ThorAxeMSA.select_best_seed(summary)
    @test larger_msa_seed.pid == 80.0

    write(summary,
        "pid,pid_order,eligible,epli,n_sequences_msa0,stockholm_path,fasta_path\n" *
        "10.0,2,true,70.0,20,pid10.sto,pid10.fa\n" *
        "80.0,1,true,70.0,20,pid80.sto,pid80.fa\n")
    ordered_seed = Iduna.ThorAxeMSA.select_best_seed(summary)
    @test ordered_seed.pid == 80.0
    Iduna.ThorAxeMSA._mark_selected_candidate!(summary, ordered_seed)
    selected_lines = read(summary, String)
    @test occursin("80.0,1,true,70.0,20,pid80.sto,pid80.fa,true",
        selected_lines)
    all_candidate_seeds = Iduna.ThorAxeMSA._select_scored_candidate_seeds(summary, 0)
    @test [seed.pid for seed in all_candidate_seeds] == [80.0, 10.0]
    sampled_seed = only(Iduna.ThorAxeMSA._select_scored_candidate_seeds(summary, 1))
    @test sampled_seed.pid == 80.0

    write(summary,
        "pid,epli,stockholm_path,fasta_path\n" *
        "80.0,70.0,pid80.sto,pid80.fa\n" *
        "10.0,70.0,pid10.sto,pid10.fa\n")
    tied_seed = Iduna.ThorAxeMSA.select_best_seed(summary)
    @test tied_seed.pid == 80.0

    write(summary,
        "pid,eligible,epli,stockholm_path,fasta_path\n" *
        "30.0,false,90.0,pid30.sto,pid30.fa\n")
    @test_throws ErrorException Iduna.ThorAxeMSA.select_best_seed(summary)
end
