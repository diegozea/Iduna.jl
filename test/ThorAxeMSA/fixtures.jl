function write_test_ensembl_bundle(root::AbstractString)
    ensembl = joinpath(root, "Ensembl")
    mkpath(ensembl)
    for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
        write(joinpath(ensembl, file), "$(file)\n")
    end
    return root
end

function write_fake_thoraxe_dir(thoraxe_root::AbstractString)
    msa_dir = joinpath(thoraxe_root, "msa")
    phylosofs_dir = joinpath(thoraxe_root, "phylosofs")
    mkpath(msa_dir)
    mkpath(phylosofs_dir)
    write(joinpath(thoraxe_root, "path_table.csv"),
        "TranscriptIDCluster,Path\nENST,start/1_0/stop\n")
    write(joinpath(thoraxe_root, "s_exon_table.csv"),
        "GeneID,Species,TranscriptIDCluster,S_exonID,S_exon_Sequence\n" *
        "ENSG,homo_sapiens,ENST,1_0,AA\n" *
        "ORTHO1,mus_musculus,ORTHO1,1_0,AB\n")
    write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
        ">ENSG\nAA\n>ORTHO1\nAB\n")
    write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0	a\n")
    write(joinpath(phylosofs_dir, "transcripts.pir"),
        ">P1;ENSG ENST a\naa\nAA*\n")
    return thoraxe_root
end
