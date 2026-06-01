@testset "ThorAxeMSA" begin
    function write_test_ensembl_bundle(root::AbstractString)
        ensembl = joinpath(root, "Ensembl")
        mkpath(ensembl)
        for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
            write(joinpath(ensembl, file), "$(file)\n")
        end
        return root
    end

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
            "pid,pid_order,eligible,median_identity,mean_identity,n_sequences_msa0,stockholm_path,fasta_path\n" *
            "30.0,1,false,90.0,90.0,20,pid30.sto,pid30.fa\n" *
            "10.0,2,true,70.0,80.0,10,pid10.sto,pid10.fa\n" *
            "80.0,3,true,60.0,90.0,30,pid80.sto,pid80.fa\n")
        seed = Iduna.ThorAxeMSA.select_best_seed(summary)
        @test seed.pid == 10.0
        @test seed.median_identity == 70.0

        write(summary,
            "pid,pid_order,eligible,median_identity,mean_identity,n_sequences_msa0,stockholm_path,fasta_path\n" *
            "10.0,1,true,70.0,80.0,10,pid10.sto,pid10.fa\n" *
            "80.0,2,true,70.0,80.0,20,pid80.sto,pid80.fa\n")
        larger_msa_seed = Iduna.ThorAxeMSA.select_best_seed(summary)
        @test larger_msa_seed.pid == 80.0

        write(summary,
            "pid,pid_order,eligible,median_identity,mean_identity,n_sequences_msa0,stockholm_path,fasta_path\n" *
            "10.0,2,true,70.0,80.0,20,pid10.sto,pid10.fa\n" *
            "80.0,1,true,70.0,80.0,20,pid80.sto,pid80.fa\n")
        ordered_seed = Iduna.ThorAxeMSA.select_best_seed(summary)
        @test ordered_seed.pid == 80.0
        Iduna.ThorAxeMSA._mark_selected_candidate!(summary, ordered_seed)
        selected_lines = read(summary, String)
        @test occursin("80.0,1,true,70.0,80.0,20,pid80.sto,pid80.fa,true",
            selected_lines)
        all_candidate_seeds = Iduna.ThorAxeMSA._select_scored_candidate_seeds(summary, 0)
        @test [seed.pid for seed in all_candidate_seeds] == [80.0, 10.0]
        sampled_seed = only(Iduna.ThorAxeMSA._select_scored_candidate_seeds(summary, 1))
        @test sampled_seed.pid == 80.0

        write(summary,
            "pid,median_identity,mean_identity,stockholm_path,fasta_path\n" *
            "80.0,70.0,80.0,pid80.sto,pid80.fa\n" *
            "10.0,70.0,80.0,pid10.sto,pid10.fa\n")
        tied_seed = Iduna.ThorAxeMSA.select_best_seed(summary)
        @test tied_seed.pid == 80.0

        write(summary,
            "pid,eligible,median_identity,mean_identity,stockholm_path,fasta_path\n" *
            "30.0,false,90.0,90.0,pid30.sto,pid30.fa\n")
        @test_throws ErrorException Iduna.ThorAxeMSA.select_best_seed(summary)
    end

    mktempdir() do tmp
        fasta = joinpath(tmp, "candidate.fasta")
        write(fasta,
            ">ORTHO1\nAAAA\n" *
            ">ENSG00000000001\nBBBB\n" *
            ">ORTHO2\nCCCC\n" *
            ">ORTHO3\nDDDD\n")
        msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
        species = ["mus_musculus", "homo_sapiens", "rattus_norvegicus", "danio_rerio"]

        rng = Iduna.ThorAxeMSA._sample_rng(UInt64(7), 1)
        indices = Iduna.ThorAxeMSA._sample_indices(4, 2, 0.5, rng)
        @test first(indices) == 2
        @test length(indices) == 3
        @test length(unique(indices)) == length(indices)
        @test !(2 in indices[2:end])

        full_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 80.0, 0)
        Iduna.ThorAxeMSA._write_candidate_sample_inputs(
            full_paths, msa, species, collect(1:4); overwrite = true)
        Iduna.ThorAxeMSA._ensure_pid_candidate_samples(tmp, 80.0, msa, species;
            sample_count = 1,
            sample_fraction = 0.5,
            sample_seed = UInt64(7),
            overwrite = true,
            gene_id = "ENSG00000000001.1",
            transcript_id = "ENST1")

        sample_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 80.0, 1)
        @test full_paths.species_file ==
              joinpath(tmp, "thoraxe_msa", "candidates", "pid_80.00", "species",
            "candidate_species_full.txt")
        @test sample_paths.species_file ==
              joinpath(tmp, "thoraxe_msa", "candidates", "pid_80.00", "species",
            "candidate_species_subset_001.txt")
        sample0 = split(chomp(read(full_paths.species_file, String)), '\n')
        sample1 = split(chomp(read(sample_paths.species_file, String)), '\n')
        @test sample0 == species
        @test first(sample1) == "homo_sapiens"
        @test length(sample1) == 3
        @test length(unique(sample1)) == length(sample1)
        @test all(item -> item in species[[1, 3, 4]], sample1[2:end])
    end

    mktempdir() do tmp
        input_dir = joinpath(tmp, "thoraxe_input")
        ensembl = joinpath(input_dir, "Ensembl")
        mkpath(ensembl)
        for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
            write(joinpath(ensembl, file), "$(file)\n")
        end
        metadata_target = Iduna.ResolvedTarget(;
            input_id = "P20963",
            input_kind = :uniprot,
            ensembl_gene_id = "ENSG",
            transcript_id = "ENST")
        metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
            input_dir, metadata_target, [10.0, 80.0];
            sample_count = 2,
            sample_fraction = 0.5,
            sample_seed = UInt64(7),
            requested_sample_seed = 7,
            effective_specieslist = "homo_sapiens",
            orthology = "1:1",
            specieslist_filter = true,
            biomart_datasets_filter = false)
        summary = joinpath(tmp, "candidate_summary.csv")
        rows = [
            (;
                gene_id = "ENSG",
                transcript_id = "ENST",
                pid = 10.0,
                pid_order = 1,
                eligible = true,
                selected = true,
                msa0_status = "ok",
                msa0_issue = missing,
                mean_identity = 70.0,
                median_identity = 70.0,
                n_samples = 2,
                n_sequences_msa0 = 4,
                pid_sample_count = metadata.pid_sample_count,
                pid_sample_fraction = metadata.pid_sample_fraction,
                pid_sample_seed = metadata.pid_sample_seed,
                pid_thresholds_key = metadata.pid_thresholds_key,
                effective_specieslist = metadata.effective_specieslist,
                orthology = metadata.orthology,
                specieslist_filter = metadata.specieslist_filter,
                biomart_datasets_filter = metadata.biomart_datasets_filter,
                transcript_query_fingerprint = metadata.transcript_query_fingerprint,
                selection_mode = metadata.selection_mode,
                fasta_path = "pid10.fasta",
                stockholm_path = "pid10.sto",
                sequence_fasta = "pid10_sequences.fasta",
                species_file = "pid10_species.txt",
                scores_path = "pid10_scores.csv"
            ),
            (;
                gene_id = "ENSG",
                transcript_id = "ENST",
                pid = 80.0,
                pid_order = 2,
                eligible = false,
                selected = false,
                msa0_status = "invalid_msa0",
                msa0_issue = "indels",
                mean_identity = missing,
                median_identity = missing,
                n_samples = 0,
                n_sequences_msa0 = 3,
                pid_sample_count = metadata.pid_sample_count,
                pid_sample_fraction = metadata.pid_sample_fraction,
                pid_sample_seed = metadata.pid_sample_seed,
                pid_thresholds_key = metadata.pid_thresholds_key,
                effective_specieslist = metadata.effective_specieslist,
                orthology = metadata.orthology,
                specieslist_filter = metadata.specieslist_filter,
                biomart_datasets_filter = metadata.biomart_datasets_filter,
                transcript_query_fingerprint = metadata.transcript_query_fingerprint,
                selection_mode = metadata.selection_mode,
                fasta_path = "pid80.fasta",
                stockholm_path = "pid80.sto",
                sequence_fasta = "pid80_sequences.fasta",
                species_file = "pid80_species.txt",
                scores_path = "pid80_scores.csv"
            )
        ]
        Iduna.ThorAxeMSA.CSV.write(summary, Iduna.ThorAxeMSA.DataFrame(rows))
        df = Iduna.ThorAxeMSA._candidate_summary_dataframe(summary)
        @test Iduna.ThorAxeMSA._candidate_summary_matches(df, metadata)
        @test Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, metadata)
        changed_pids = merge(metadata, (; pid_thresholds_key = "80.0"))
        @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_pids)
        @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_pids)
        changed_fraction = merge(metadata, (; pid_sample_fraction = 0.75))
        @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_fraction)
        @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_fraction)
        changed_gene = merge(metadata, (; gene_id = "ENSG_DIFFERENT"))
        @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_gene)
        @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_gene)
        changed_transcript = merge(metadata, (; transcript_id = "ENST_DIFFERENT"))
        @test !Iduna.ThorAxeMSA._candidate_summary_matches(df, changed_transcript)
        @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, changed_transcript)
        @test !Iduna.ThorAxeMSA._has_matching_candidate_summary(
            joinpath(tmp, "missing_summary.csv"), metadata)
        random_requested_seed = merge(metadata,
            (; pid_sample_seed = UInt64(99), requested_pid_sample_seed = nothing))
        @test Iduna.ThorAxeMSA._candidate_summary_matches(df, random_requested_seed)
        @test Iduna.ThorAxeMSA._has_matching_candidate_summary(summary, random_requested_seed)
        @test Iduna.ThorAxeMSA._thoraxe_msa_identity(metadata, tmp) !=
              Iduna.ThorAxeMSA._thoraxe_msa_identity(
            merge(metadata,
                (; pid_sample_seed = UInt64(99), requested_pid_sample_seed = 99)),
            tmp)
        @test Iduna.ThorAxeMSA._thoraxe_msa_identity(random_requested_seed, tmp) ==
              Iduna.ThorAxeMSA._thoraxe_msa_identity(
            merge(random_requested_seed, (; pid_sample_seed = UInt64(7))), tmp)
        random_identity = Iduna.ThorAxeMSA._thoraxe_msa_identity(
            random_requested_seed, tmp)
        Iduna.ThorAxeMSA._write_thoraxe_msa_state(
            tmp, summary, Iduna.SeedSelection[], :done, random_identity; action = :run)
        random_stage_cache = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(
            tmp, summary, random_requested_seed,
            Iduna.ThorAxeMSA._thoraxe_msa_identity(
                merge(random_requested_seed, (; pid_sample_seed = UInt64(7))), tmp);
            overwrite = false)
        @test random_stage_cache.cache.reusable === true

        no_species_metadata = merge(metadata, (; effective_specieslist = nothing))
        fasta = joinpath(tmp, "candidate_summary_candidate.fasta")
        write(fasta,
            ">ENSG\nAAAA\n" *
            ">ORTHO1\nBBBB\n")
        candidate_msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
        candidate = (;
            msa = candidate_msa,
            fasta_path = "pid10.fasta",
            stockholm_path = "pid10.sto",
            sequence_fasta = "pid10_sequences.fasta",
            species_file = "pid10_species.txt",
            workdir = tmp)
        validation = (; eligible = true, status = "ok", issue = missing)
        no_species_rows = [
            Iduna.ThorAxeMSA._candidate_summary_row(
                metadata_target, candidate, 10.0, 1, validation;
                sample_count = no_species_metadata.pid_sample_count,
                sample_fraction = no_species_metadata.pid_sample_fraction,
                sample_seed = no_species_metadata.pid_sample_seed,
                metadata = no_species_metadata),
            Iduna.ThorAxeMSA._candidate_summary_row(
                metadata_target, candidate, 80.0, 2, validation;
                sample_count = no_species_metadata.pid_sample_count,
                sample_fraction = no_species_metadata.pid_sample_fraction,
                sample_seed = no_species_metadata.pid_sample_seed,
                metadata = no_species_metadata)
        ]
        no_species_summary = joinpath(tmp, "candidate_summary_no_species.csv")
        Iduna.ThorAxeMSA._summarize_candidate_scores(
            no_species_rows, no_species_summary)
        no_species_df = Iduna.ThorAxeMSA._candidate_summary_dataframe(no_species_summary)
        @test all(ismissing, no_species_df.effective_specieslist)
        @test Iduna.ThorAxeMSA._candidate_summary_matches(
            no_species_df, no_species_metadata)

        decimal_pid_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
            input_dir, metadata_target, [12.34];
            sample_count = 2,
            sample_fraction = 0.5,
            sample_seed = UInt64(7),
            requested_sample_seed = 7,
            effective_specieslist = "homo_sapiens",
            orthology = "1:1",
            specieslist_filter = true,
            biomart_datasets_filter = false)
        @test decimal_pid_metadata.pid_thresholds_key == "12.34"
        decimal_pid_rows = [
            Iduna.ThorAxeMSA._candidate_summary_row(
            metadata_target, candidate, 12.34, 1, validation;
            mean_identity = 70.0,
            median_identity = 70.0,
            n_samples = 2,
            sample_count = decimal_pid_metadata.pid_sample_count,
            sample_fraction = decimal_pid_metadata.pid_sample_fraction,
            sample_seed = decimal_pid_metadata.pid_sample_seed,
            metadata = decimal_pid_metadata)
        ]
        decimal_pid_summary = joinpath(tmp, "candidate_summary_decimal_pid.csv")
        Iduna.ThorAxeMSA._summarize_candidate_scores(
            decimal_pid_rows, decimal_pid_summary)
        decimal_pid_df = Iduna.ThorAxeMSA._candidate_summary_dataframe(
            decimal_pid_summary)
        @test Iduna.ThorAxeMSA._candidate_summary_matches(
            decimal_pid_df, decimal_pid_metadata)

        selected = Iduna.ThorAxeMSA._selected_candidate_seeds(summary)
        @test length(selected) == 1
        @test only(selected).pid == 10.0
        eligible = Iduna.ThorAxeMSA._eligible_candidate_seeds(summary)
        @test length(eligible) == 1
        @test only(eligible).pid == 10.0

        warning_summary = joinpath(tmp, "candidate_summary_warnings.csv")
        warning_rows = [
            Iduna.ThorAxeMSA._candidate_summary_row(
                metadata_target, candidate, 20.0, 1,
                (; eligible = false, status = "invalid_msa0", issue = "bad indel");
                sample_count = metadata.pid_sample_count,
                sample_fraction = metadata.pid_sample_fraction,
                sample_seed = metadata.pid_sample_seed,
                metadata),
            Iduna.ThorAxeMSA._candidate_summary_row(
                metadata_target, candidate, 30.0, 2,
                (; eligible = true, status = "warning", issue = "one substitution");
                sample_count = metadata.pid_sample_count,
                sample_fraction = metadata.pid_sample_fraction,
                sample_seed = metadata.pid_sample_seed,
                metadata)
        ]
        Iduna.ThorAxeMSA._summarize_candidate_scores(warning_rows, warning_summary)
        warning_text = Iduna.ThorAxeMSA._candidate_summary_warnings(warning_summary)
        @test any(w -> occursin("excluded", w) && occursin("bad indel", w),
            warning_text)
        @test any(
            w -> occursin("candidate retained", w) &&
                 occursin("one substitution", w), warning_text)
    end

    mktempdir() do tmp
        stdout_log = joinpath(tmp, "stdout.log")
        stderr_log = joinpath(tmp, "stderr.log")
        Iduna.ThorAxeMSA._run_logged_command(
            `$(Base.julia_cmd()) --startup-file=no -e "println(\"ok\")"`,
            stdout_log,
            stderr_log)
        @test read(stdout_log, String) == "ok\n"
        @test isempty(read(stderr_log, String))

        @test_throws ErrorException Iduna.ThorAxeMSA._run_logged_command(
            `sh -c "exit 2"`,
            joinpath(tmp, "failed_stdout.log"),
            joinpath(tmp, "failed_stderr.log"))
    end

    mktempdir() do tmp
        source = joinpath(tmp, "cached_input")
        ensembl = joinpath(source, "Ensembl")
        mkpath(ensembl)
        for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
            write(joinpath(ensembl, file), "x\n")
        end
        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1")
        workdir = joinpath(tmp, "work")
        copied = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
            cached_input_dir = source,
            overwrite = true)
        @test copied == joinpath(workdir, "thoraxe_input")
        @test isfile(joinpath(copied, "Ensembl", "sequences.fasta"))
        @test isdir(source)
        reused = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
            cached_input_dir = source,
            overwrite = false)
        @test reused == copied
        write(joinpath(source, "Ensembl", "sequences.fasta"), "changed\n")
        recopied = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
            cached_input_dir = source,
            overwrite = false)
        @test recopied == copied
        @test read(joinpath(recopied, "Ensembl", "sequences.fasta"), String) == "changed\n"

        @test_throws ErrorException Iduna.ThorAxeMSA._ensure_transcript_query(
            target, joinpath(tmp, "bad_work");
            cached_input_dir = joinpath(tmp, "missing"),
            overwrite = true)

        invalid_metadata_source = joinpath(tmp, "invalid_metadata_source")
        write_test_ensembl_bundle(invalid_metadata_source)
        mkpath(joinpath(invalid_metadata_source,
            Iduna.ThorAxeMSA._TRANSCRIPT_QUERY_METADATA_FILE))
        failed_copy_workdir = joinpath(tmp, "failed_copy_work")
        @test_throws SystemError Iduna.ThorAxeMSA._ensure_transcript_query(
            target, failed_copy_workdir;
            cached_input_dir = invalid_metadata_source,
            overwrite = true)
        failed_input_state = Iduna.Utils._read_stage_state(
            Iduna.ThorAxeMSA._thoraxe_input_stage_dir(failed_copy_workdir))
        @test failed_input_state["status"] == "failed"
        @test failed_input_state["exception"]["type"] == "SystemError"

        direct_workdir = joinpath(tmp, "direct_work")
        direct_input = Iduna.ThorAxeMSA._thoraxe_input_dir(direct_workdir)
        write_test_ensembl_bundle(direct_input)
        direct_metadata = Iduna.ThorAxeMSA._expected_transcript_query_metadata(target;
            specieslist = "homo_sapiens",
            orthology = "1:1",
            source_kind = "transcript_query")
        Iduna.ThorAxeMSA._write_transcript_query_metadata!(
            direct_input, direct_metadata)
        @test Iduna.ThorAxeMSA._ensure_transcript_query(target, direct_workdir;
            specieslist = "homo_sapiens",
            orthology = "1:1") == direct_input

        queried_workdir = joinpath(tmp, "queried_work")
        mkpath(queried_workdir)
        query_calls = Ref(0)
        query_runner = command -> begin
            query_calls[] += 1
            write_test_ensembl_bundle(joinpath(pwd(), "ENSG00000000001"))
            nothing
        end
        rebuilt_input = Iduna.ThorAxeMSA._ensure_transcript_query(
            target, queried_workdir;
            specieslist = "homo_sapiens",
            orthology = "1:1",
            transcript_query_runner = query_runner,
            sleep_fn = seconds -> nothing)
        @test rebuilt_input == Iduna.ThorAxeMSA._thoraxe_input_dir(queried_workdir)
        @test query_calls[] == 1

        default_runner_workdir = joinpath(tmp, "default_runner_work")
        mkpath(default_runner_workdir)
        factory_calls = Ref(0)
        runner_calls = Ref(0)
        captured_logs = Tuple{String, String}[]
        default_runner_factory = (stdout_log,
            stderr_log) -> begin
            factory_calls[] += 1
            push!(captured_logs, (stdout_log, stderr_log))
            return command -> begin
                runner_calls[] += 1
                write_test_ensembl_bundle(joinpath(pwd(), "ENSG00000000001"))
                nothing
            end
        end
        default_runner_input = Iduna.ThorAxeMSA._ensure_transcript_query(
            target, default_runner_workdir;
            specieslist = "homo_sapiens",
            orthology = "1:1",
            thoraxe_runner_factory = default_runner_factory,
            sleep_fn = seconds -> nothing)
        @test default_runner_input ==
              Iduna.ThorAxeMSA._thoraxe_input_dir(default_runner_workdir)
        @test factory_calls[] == 1
        @test runner_calls[] == 1
        @test only(captured_logs) == (
            joinpath(default_runner_workdir, "logs", "thoraxe",
                "transcript_query_stdout.log"),
            joinpath(default_runner_workdir, "logs", "thoraxe",
                "transcript_query_stderr.log"))
    end

    @testset "metadata and cached branch helpers" begin
        mktempdir() do tmp
            metadata_path = joinpath(tmp, "metadata.json")
            expected = (gene_id = "ENSG", specieslist = nothing)
            @test !Iduna.ThorAxeMSA._metadata_matches(metadata_path, expected)

            write(metadata_path, """{"gene_id":"ENSG","specieslist":null}""")
            @test Iduna.ThorAxeMSA._metadata_matches(metadata_path, expected)
            @test !Iduna.ThorAxeMSA._metadata_matches(
                metadata_path, (; gene_id = "OTHER", specieslist = nothing))
            @test !Iduna.ThorAxeMSA._metadata_matches(
                metadata_path, (; gene_id = "ENSG", transcript_id = "ENST"))

            write(metadata_path, "{not json")
            @test !Iduna.ThorAxeMSA._metadata_matches(metadata_path, expected)
            @test Iduna.ThorAxeMSA._retry_wait_seconds(6) == 30.0
        end

        mktempdir() do tmp
            fasta = joinpath(tmp, "candidate.fasta")
            write(fasta, ">ENSG\nA-C\n>ORTHO1\nABC\n")
            msa = Iduna.ThorAxeMSA.read_file(fasta, Iduna.ThorAxeMSA.FASTA)
            Iduna.Utils.set_s_exon_annotations!(msa, "000", ['0' => "1_0"])
            paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 60.0, 0)
            Iduna.ThorAxeMSA._write_candidate_sample_inputs(
                paths, msa, ["homo_sapiens", "mus_musculus"], [1, 2])
            Iduna.ThorAxeMSA.write_file(paths.fasta_path, msa, Iduna.ThorAxeMSA.FASTA)
            Iduna.ThorAxeMSA.write_file(
                paths.stockholm_path, msa, Iduna.ThorAxeMSA.Stockholm)
            written_species = read(paths.species_file, String)
            @test Iduna.ThorAxeMSA._write_species_file(
                paths.species_file, ["danio_rerio"]) == paths.species_file
            @test read(paths.species_file, String) == written_species
            @test Iduna.ThorAxeMSA._write_candidate_sample_inputs(
                paths, msa, ["wrong_length"], [1]) ==
                  (paths.sequence_fasta, paths.species_file)
            @test_throws ErrorException Iduna.ThorAxeMSA._write_candidate_sample_inputs(
                Iduna.ThorAxeMSA._pid_sample_paths(tmp, 70.0, 0),
                msa, ["wrong_length"], [1])

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
                write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\ta\n")
                write(joinpath(phylosofs_dir, "transcripts.pir"),
                    ">P1;ENSG ENST a\naa\nAA*\n")
                return thoraxe_root
            end

            thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 60.0, 0)
            mkpath(thoraxe_dir)
            write(joinpath(thoraxe_dir, "path_table.csv"), "TranscriptIDCluster,Path\n")
            @test !Iduna.ThorAxeMSA._has_phylosofs_outputs(thoraxe_dir)
            phylosofs_dir = joinpath(thoraxe_dir, "phylosofs")
            mkpath(phylosofs_dir)
            write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\t0\n")
            write(joinpath(phylosofs_dir, "transcripts.pir"), ">P1;ENST\n0\nA*\n")
            @test Iduna.ThorAxeMSA._has_phylosofs_outputs(thoraxe_dir)
            @test Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
                Iduna.ResolvedTarget(;
                    input_id = "ENST",
                    input_kind = :ensembl_transcript,
                    ensembl_gene_id = "ENSG",
                    transcript_id = "ENST"),
                joinpath(tmp, "input"), tmp, 60.0, nothing, 0;
                keep_thoraxe_dir = true) ==
                  (paths.fasta_path, paths.stockholm_path, thoraxe_dir)

            kept_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 90.0, 0)
            kept_thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 90.0, 0)
            mkpath(dirname(kept_paths.fasta_path))
            write_fake_thoraxe_dir(kept_thoraxe_dir)
            @test Iduna.ThorAxeMSA._run_kept_thoraxe_pid_msa!(
                Iduna.ResolvedTarget(;
                    input_id = "ENST",
                    input_kind = :ensembl_transcript,
                    ensembl_gene_id = "ENSG",
                    transcript_id = "ENST"),
                joinpath(tmp, "input"), tmp, kept_paths, kept_thoraxe_dir,
                joinpath(kept_thoraxe_dir, "path_table.csv"), 90.0, nothing, 0,
                command -> nothing, false) ==
                  (kept_paths.fasta_path, kept_paths.stockholm_path, kept_thoraxe_dir)
            @test isfile(kept_paths.stockholm_path)

            fake_calls = Ref(0)
            fake_thoraxe = (input_dir, run_root; identity, specieslist, phylosofs,
                runner) -> begin
                fake_calls[] += 1
                @test identity == 50.0
                @test phylosofs === true
                write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
                nothing
            end
            temp_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 50.0, 0)
            @test Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
                Iduna.ResolvedTarget(;
                    input_id = "ENST",
                    input_kind = :ensembl_transcript,
                    ensembl_gene_id = "ENSG",
                    transcript_id = "ENST"),
                joinpath(tmp, "input"), tmp, 50.0, nothing, 0;
                thoraxe_fn = fake_thoraxe) ==
                  (temp_paths.fasta_path, temp_paths.stockholm_path,
                Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 50.0, 0))
            @test fake_calls[] == 1
            @test isfile(temp_paths.stockholm_path)

            generated_pid = 55.0
            generated_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, generated_pid, 0)
            generated_thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(
                tmp, generated_pid, 0)
            write_fake_thoraxe_dir(generated_thoraxe_dir)
            generated_candidate = Iduna.ThorAxeMSA._generate_pid_candidate(
                Iduna.ResolvedTarget(;
                    input_id = "ENST",
                    input_kind = :ensembl_transcript,
                    ensembl_gene_id = "ENSG",
                    transcript_id = "ENST"),
                joinpath(tmp, "input"), tmp, generated_pid, nothing)
            @test generated_candidate.fasta_path == generated_paths.fasta_path
            @test generated_candidate.stockholm_path == generated_paths.stockholm_path
            @test generated_candidate.thoraxe_dir == generated_thoraxe_dir
            @test generated_candidate.s_exon_blocks_tsv == generated_paths.s_exon_blocks_tsv
            @test isfile(generated_candidate.stockholm_path)
            @test isfile(generated_candidate.s_exon_blocks_tsv)
            generated_sequence_fasta = read(generated_candidate.sequence_fasta, String)
            @test occursin(">ENSG\nAA", generated_sequence_fasta)
            @test occursin(">ORTHO1\nAX", generated_sequence_fasta)
            @test read(generated_candidate.species_file, String) ==
                  "homo_sapiens\nmus_musculus\n"

            scoring_target = Iduna.ResolvedTarget(;
                input_id = "ENST",
                input_kind = :ensembl_transcript,
                ensembl_gene_id = "ENSG",
                transcript_id = "ENST")
            scoring_input = write_test_ensembl_bundle(joinpath(tmp, "scoring_input"))
            zero_sample_pid = 52.0
            write_fake_thoraxe_dir(
                Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, zero_sample_pid, 0))
            zero_sample_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
                scoring_input, scoring_target, [zero_sample_pid];
                sample_count = 0,
                sample_fraction = 1.0,
                sample_seed = UInt64(11),
                requested_sample_seed = 11,
                effective_specieslist = nothing,
                orthology = "1:1",
                specieslist_filter = false,
                biomart_datasets_filter = false)
            zero_sample_row = Iduna.ThorAxeMSA._score_pid_candidate(
                scoring_target, scoring_input, tmp, zero_sample_pid, 1, nothing;
                sample_count = 0,
                sample_fraction = 1.0,
                sample_seed = UInt64(11),
                metadata = zero_sample_metadata)
            @test zero_sample_row.pid == zero_sample_pid
            @test zero_sample_row.eligible
            @test zero_sample_row.selection_mode == "all_candidates"
            @test zero_sample_row.n_samples == 0

            invalid_uniprot = joinpath(tmp, "invalid_uniprot.fasta")
            write(invalid_uniprot, ">P0\nA\n")
            invalid_target = Iduna.ResolvedTarget(;
                input_id = "P0",
                input_kind = :uniprot,
                uniprot_id = "P0",
                ensembl_gene_id = "ENSG",
                transcript_id = "ENST",
                uniprot_sequence_path = invalid_uniprot)
            invalid_pid = 52.5
            write_fake_thoraxe_dir(
                Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, invalid_pid, 0))
            invalid_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
                scoring_input, invalid_target, [invalid_pid];
                sample_count = 0,
                sample_fraction = 1.0,
                sample_seed = UInt64(13),
                requested_sample_seed = 13,
                effective_specieslist = nothing,
                orthology = "1:1",
                specieslist_filter = false,
                biomart_datasets_filter = false)
            invalid_row = @test_logs (:info, r"Skipping ineligible ThorAxe PID candidate") match_mode=:any Iduna.ThorAxeMSA._score_pid_candidate(
                invalid_target, scoring_input, tmp, invalid_pid, 1, nothing;
                sample_count = 0,
                sample_fraction = 1.0,
                sample_seed = UInt64(13),
                metadata = invalid_metadata)
            @test invalid_row.pid == invalid_pid
            @test !invalid_row.eligible
            @test invalid_row.msa0_status == "invalid_msa0"
            @test occursin("indels", invalid_row.msa0_issue)

            sampled_pid = 53.0
            sampled_calls = Ref(0)
            sampled_thoraxe = (input_dir, run_root; identity, specieslist, phylosofs,
                runner) -> begin
                sampled_calls[] += 1
                write_fake_thoraxe_dir(joinpath(run_root, "thoraxe"))
                nothing
            end
            sampled_identity = (reference_fasta, sample_fasta; logs_dir,
                label) -> begin
                @test isfile(reference_fasta)
                @test isfile(sample_fasta)
                @test label == "sample1"
                mkpath(logs_dir)
                write(joinpath(logs_dir, "$(label)_hhalign.out"), "fake\n")
                88.0
            end
            sampled_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
                scoring_input, scoring_target, [sampled_pid];
                sample_count = 1,
                sample_fraction = 1.0,
                sample_seed = UInt64(12),
                requested_sample_seed = 12,
                effective_specieslist = nothing,
                orthology = "1:1",
                specieslist_filter = false,
                biomart_datasets_filter = false)
            sampled_row = Iduna.ThorAxeMSA._score_pid_candidate(
                scoring_target, scoring_input, tmp, sampled_pid, 1, nothing;
                sample_count = 1,
                sample_fraction = 1.0,
                sample_seed = UInt64(12),
                metadata = sampled_metadata,
                thoraxe_fn = sampled_thoraxe,
                identity_fn = sampled_identity)
            @test sampled_calls[] == 2
            @test sampled_row.mean_identity == 88.0
            @test sampled_row.median_identity == 88.0
            @test sampled_row.n_samples == 1
            @test isfile(Iduna.ThorAxeMSA._pid_scores_path(tmp, sampled_pid))

            scored_rows = Iduna.ThorAxeMSA._score_pid_candidates(
                scoring_target, scoring_input, tmp, [62.0, 63.0], nothing,
                zero_sample_metadata;
                pid_sample_count = 0,
                pid_sample_fraction = 1.0,
                sample_seed = UInt64(13),
                overwrite = false,
                thoraxe_fn = sampled_thoraxe,
                identity_fn = sampled_identity)
            @test [row.pid for row in scored_rows] == [62.0, 63.0]
            @test all(row -> row.selection_mode == "all_candidates", scored_rows)

            stale_paths = Iduna.ThorAxeMSA._pid_sample_paths(tmp, 80.0, 0)
            mkpath(dirname(stale_paths.fasta_path))
            Iduna.ThorAxeMSA.write_file(stale_paths.fasta_path, msa, Iduna.ThorAxeMSA.FASTA)
            Iduna.ThorAxeMSA.write_file(
                stale_paths.stockholm_path, msa, Iduna.ThorAxeMSA.Stockholm)
            stale_thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(tmp, 80.0, 0)
            mkpath(stale_thoraxe_dir)
            write(joinpath(stale_thoraxe_dir, "path_table.csv"), "TranscriptIDCluster,Path\n")
            @test_throws Exception Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
                Iduna.ResolvedTarget(;
                    input_id = "ENST",
                    input_kind = :ensembl_transcript,
                    ensembl_gene_id = "ENSG",
                    transcript_id = "ENST"),
                joinpath(tmp, "missing_input"), tmp, 80.0, nothing, 0;
                keep_thoraxe_dir = true)
            @test_throws Exception Iduna.ThorAxeMSA._run_thoraxe_pid_msa(
                Iduna.ResolvedTarget(;
                    input_id = "ENST",
                    input_kind = :ensembl_transcript,
                    ensembl_gene_id = "ENSG",
                    transcript_id = "ENST"),
                joinpath(tmp, "missing_input"), tmp, 60.0, nothing, 0;
                overwrite = true,
                keep_thoraxe_dir = true)

            positions,
            codes = Iduna.ThorAxeMSA._get_codes(
                "Probab=99\nquery 1 A-C 2\n        | |\nConfidence\n")
            @test positions == [1, 0, 2]
            @test codes == ['|', ' ', '|']
            @test Iduna.ThorAxeMSA._identity_from_codes(positions, codes) == 100.0
            @test_throws ErrorException Iduna.ThorAxeMSA._reference_index(
                msa, "MISSING_GENE", "MISSING_TRANSCRIPT")

            reference_fasta = joinpath(tmp, "identity_reference.fasta")
            sample_fasta = joinpath(tmp, "identity_sample.fasta")
            write(reference_fasta, ">ref\nAA\n")
            write(sample_fasta, ">sample\nAA\n")
            identity_logs = joinpath(tmp, "identity_logs")
            @test Iduna.ThorAxeMSA.compute_identity_against_reference(
                reference_fasta, sample_fasta; logs_dir = identity_logs,
                label = "identical") == 100.0
            @test isfile(joinpath(identity_logs, "identical_hhalign.out"))
        end
    end

    @testset "orthology specieslist filter helpers" begin
        @test Iduna.ThorAxeMSA._normalize_species_name(nothing) === nothing
        @test Iduna.ThorAxeMSA._prepend_query_species(["mus_musculus"], nothing) ==
              ["mus_musculus"]
        @test Iduna.ThorAxeMSA._orthology_relationships("1:1") ==
              ["ortholog_one2one"]
        @test Iduna.ThorAxeMSA._orthology_relationships("1:n") ==
              ["ortholog_one2one", "ortholog_one2many"]
        @test Iduna.ThorAxeMSA._orthology_relationships("m:n") ==
              ["ortholog_one2one", "ortholog_one2many", "ortholog_many2many"]
        @test_throws ErrorException Iduna.ThorAxeMSA._orthology_relationships("all")

        homologies = [
            Dict(
                "type" => "ortholog_one2one",
                "target" => Dict("species" => "mus_musculus")),
            Dict(
                "type" => "ortholog_one2many",
                "target" => Dict("species" => "danio_rerio")),
            Dict(
                "type" => "ortholog_many2many",
                "target" => Dict("species" => "xenopus_tropicalis"))
        ]
        homology_data = Dict("data" => [Dict("homologies" => homologies)])
        @test Iduna.ThorAxeMSA._homology_species(homology_data, "1:1") ==
              ["mus_musculus"]
        @test Iduna.ThorAxeMSA._homology_species(homology_data, "1:n") ==
              ["mus_musculus", "danio_rerio"]
        @test Iduna.ThorAxeMSA._homology_species(homology_data, "m:n") ==
              ["mus_musculus", "danio_rerio", "xenopus_tropicalis"]

        unknown_species_target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1")
        @test_throws ErrorException Iduna.ThorAxeMSA._fetch_ortholog_species(
            unknown_species_target, "1:1")

        species_target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1",
            species = "Homo sapiens")
        fetched_species = Iduna.ThorAxeMSA._fetch_ortholog_species(
            species_target, "1:n";
            homology_data_fetcher = (
                species, gene_id) -> begin
                @test species == "homo_sapiens"
                @test gene_id == "ENSG00000000001.1"
                homology_data
            end)
        @test fetched_species == ["mus_musculus", "danio_rerio"]
    end

    @testset "species list parsing and filtering" begin
        mktempdir() do tmp
            species_file = joinpath(tmp, "species.txt")
            write(species_file, "Mus musculus\n\nDanio rerio\n")
            @test Iduna.ThorAxeMSA._parse_specieslist(nothing) === nothing
            @test Iduna.ThorAxeMSA._parse_specieslist(
                "Homo sapiens,mus_musculus, Mus musculus ") ==
                  ["homo_sapiens", "mus_musculus"]
            @test Iduna.ThorAxeMSA._parse_specieslist(species_file) ==
                  ["mus_musculus", "danio_rerio"]
            @test Iduna.ThorAxeMSA._parse_specieslist("Canis lupus") ==
                  ["canis_lupus"]
        end

        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1",
            species = "Homo sapiens")
        fetcher = (target, orthology) -> ["mus_musculus", "danio_rerio"]

        default = Iduna.ThorAxeMSA._resolve_effective_specieslist(
            target, nothing, "1:1"; homology_species_fetcher = fetcher)
        @test default.specieslist == "homo_sapiens,mus_musculus,danio_rerio"
        @test isempty(default.warnings)

        filtered = Iduna.ThorAxeMSA._resolve_effective_specieslist(
            target, "Mus musculus,Canis lupus", "1:1";
            homology_species_fetcher = fetcher)
        @test filtered.specieslist == "homo_sapiens,mus_musculus"
        @test length(filtered.warnings) == 1
        @test occursin("canis_lupus", only(filtered.warnings))

        @test_throws ErrorException Iduna.ThorAxeMSA._resolve_effective_specieslist(
            target, "Canis lupus", "1:1"; homology_species_fetcher = fetcher)
        @test_throws ErrorException Iduna.ThorAxeMSA._resolve_effective_specieslist(
            target, nothing, "1:1";
            homology_species_fetcher = (target, orthology) -> String[])

        fallback = Iduna.ThorAxeMSA._resolve_effective_specieslist(
            target, "Canis lupus", "1:1";
            homology_species_fetcher = (target, orthology) -> error("temporary failure"))
        @test fallback.specieslist == "Canis lupus"
        @test length(fallback.warnings) == 1
        @test occursin("Ensembl specieslist filter failed", only(fallback.warnings))

        resolved_filters = Iduna.ThorAxeMSA._resolve_thoraxe_species_filters(
            target, "Mus musculus", "1:1", nothing, true, true;
            specieslist_resolver = (target,
                specieslist,
                orthology) -> (
                specieslist = "homo_sapiens,mus_musculus",
                warnings = ["species filter warning"]),
            biomart_resolver = (target,
                specieslist) -> (
                specieslist = specieslist,
                warnings = ["biomart filter warning"]))
        @test resolved_filters.effective_specieslist == "homo_sapiens,mus_musculus"
        @test resolved_filters.species_filter.warnings == ["species filter warning"]
        @test resolved_filters.biomart_filter.warnings == ["biomart filter warning"]
    end

    @testset "BioMart datasets filter" begin
        biomart_text = """

TableSet\thsapiens_gene_ensembl\tHuman genes (GRCh38.p14)\t1
TableSet\tmmusculus_gene_ensembl\tMouse genes (GRCm39)\t1
TableSet\tptroglodytes_gene_ensembl\tChimpanzee genes (Pan_tro_3.0)\t1
"""
        datasets = Iduna.ThorAxeMSA._parse_biomart_gene_datasets(biomart_text)
        @test "hsapiens_gene_ensembl" in datasets
        @test "mmusculus_gene_ensembl" in datasets
        @test "ptroglodytes_gene_ensembl" in datasets
        @test Iduna.ThorAxeMSA._biomart_gene_dataset_for_species("Homo sapiens") ==
              "hsapiens_gene_ensembl"
        @test Iduna.ThorAxeMSA._biomart_gene_dataset_for_species("mus_musculus") ==
              "mmusculus_gene_ensembl"

        target = Iduna.ResolvedTarget(;
            input_id = "ENST00000000001.1",
            input_kind = :ensembl_transcript,
            ensembl_gene_id = "ENSG00000000001.1",
            transcript_id = "ENST00000000001.1",
            species = "Homo sapiens")
        loader = () -> (
            datasets = Set(["hsapiens_gene_ensembl", "mmusculus_gene_ensembl"]),
            warnings = String[])
        filtered = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
            target, "Homo sapiens,Mus musculus,Canis lupus";
            dataset_loader = loader)
        @test filtered.specieslist == "homo_sapiens,mus_musculus"
        @test length(filtered.warnings) == 1
        @test occursin("canis_lupus", only(filtered.warnings))

        missing_query = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
            target, "Homo sapiens,Mus musculus";
            dataset_loader = () -> (datasets = Set(["mmusculus_gene_ensembl"]),
                warnings = String[]))
        @test missing_query.specieslist == "homo_sapiens,mus_musculus"
        @test any(w -> occursin("query species homo_sapiens", w),
            missing_query.warnings)

        unchecked = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
            target, "homo_sapiens,canis_lupus_familiaris";
            dataset_loader = loader)
        @test unchecked.specieslist == "homo_sapiens,canis_lupus_familiaris"
        @test any(w -> occursin("possible species aliases", w), unchecked.warnings)

        fallback = Iduna.ThorAxeMSA._resolve_biomart_datasets_specieslist(
            target, "Homo sapiens,Mus musculus";
            dataset_loader = () -> (datasets = nothing,
                warnings = ["BioMart datasets metadata refresh failed"]))
        @test fallback.specieslist == "Homo sapiens,Mus musculus"
        @test only(fallback.warnings) == "BioMart datasets metadata refresh failed"
    end

    @testset "BioMart datasets dated cache" begin
        response(status::Integer,
            body::AbstractString = "") = Iduna.ThorAxeMSA.HTTP.Response(status, Vector{UInt8}(body))

        @test isdir(Iduna.ThorAxeMSA._biomart_cache_dir())

        attempts = Ref(0)
        biomart_text = "TableSet\thsapiens_gene_ensembl\tHuman genes\t1\n"
        text = Iduna.ThorAxeMSA._fetch_biomart_datasets_text(;
            retries = 3,
            http_get = (url;
                kwargs...) -> begin
                attempts[] += 1
                attempts[] == 1 ? response(503) : response(200, biomart_text)
            end)
        @test text == biomart_text
        @test attempts[] == 2

        attempts[] = 0
        @test_throws ErrorException Iduna.ThorAxeMSA._fetch_biomart_datasets_text(;
            retries = 3,
            http_get = (url; kwargs...) -> begin
                attempts[] += 1
                response(404)
            end)
        @test attempts[] == 1

        mktempdir() do tmp
            biomart_text = "TableSet\thsapiens_gene_ensembl\tHuman genes\t1\n"
            fetches = Ref(0)
            fetcher = () -> begin
                fetches[] += 1
                biomart_text
            end
            loaded = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
                cache_dir = tmp,
                today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 6),
                fetcher)
            @test fetches[] == 1
            @test "hsapiens_gene_ensembl" in loaded.datasets
            @test isempty(loaded.warnings)

            same_day = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
                cache_dir = tmp,
                today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 6),
                fetcher = () -> error("should not refetch"))
            @test fetches[] == 1
            @test "hsapiens_gene_ensembl" in same_day.datasets

            stale = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
                cache_dir = tmp,
                today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 7),
                fetcher = () -> error("temporary failure"))
            @test "hsapiens_gene_ensembl" in stale.datasets
            @test length(stale.warnings) == 1
            @test occursin("stale cache from 2026-05-06", only(stale.warnings))
        end

        mktempdir() do tmp
            failed = Iduna.ThorAxeMSA._load_biomart_gene_datasets(;
                cache_dir = tmp,
                today = Iduna.ThorAxeMSA.Dates.Date(2026, 5, 6),
                fetcher = () -> error("temporary failure"))
            @test failed.datasets === nothing
            @test length(failed.warnings) == 1
            @test occursin("using the unfiltered specieslist", only(failed.warnings))
        end

        mktempdir() do tmp
            metadata_path = joinpath(tmp, Iduna.ThorAxeMSA._BIOMART_DATASETS_METADATA_FILE)
            @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) === nothing
            write(metadata_path, """{"download_date":"2026-05-06"}""")
            @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) == "2026-05-06"
            write(metadata_path, """{"download_time":"2026-05-06T12:00:00"}""")
            @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) === nothing
            write(metadata_path, "{not json")
            @test Iduna.ThorAxeMSA._read_biomart_cache_date(metadata_path) === nothing
        end
    end

    @testset "BioMart transcript_query warnings" begin
        mktempdir() do tmp
            input_dir = joinpath(tmp, "thoraxe_input")
            ensembl_dir = joinpath(input_dir, "Ensembl")
            logs_dir = joinpath(tmp, "logs", "thoraxe")
            mkpath(ensembl_dir)
            mkpath(logs_dir)
            write(joinpath(ensembl_dir, "errors.csv"),
                "Species,GeneID\nmus_spretus,ENSMSPG00010016579\n")
            write(joinpath(logs_dir, "transcript_query_stderr.log"),
                """
                transcript_query.py:304: UserWarning: It can not found nomascus_leucogenys in biomart (tried: ['nleucogenys_gene_ensembl', 'nleucogenys_eg_gene']).
                Last response:
                Query ERROR: caught BioMart::Exception::Usage: Dataset nleucogenys_eg_gene NOT FOUND
                  warnings.warn(...)
                transcript_query.py:728: UserWarning: Download failed for ENSPTRG00000006744 in pan_troglodytes!
                  warnings.warn(...)
                transcript_query.py:728: UserWarning: Download failed for ENSMSPG00010016579 in mus_spretus!
                  warnings.warn(...)
                """)
            warnings = Iduna.ThorAxeMSA._biomart_transcript_query_warnings(
                input_dir, logs_dir)
            @test length(warnings) == 1
            @test occursin("mus_spretus", only(warnings))
            @test occursin("nomascus_leucogenys", only(warnings))
            @test occursin("pan_troglodytes", only(warnings))
            @test occursin("errors.csv", only(warnings))
            @test occursin("transcript_query_stderr.log", only(warnings))

            malformed_errors = joinpath(input_dir, "Ensembl", "malformed_errors.csv")
            write(malformed_errors, "Other\nmus_spretus\n")
            @test isempty(Iduna.ThorAxeMSA._species_from_biomart_errors_file(
                malformed_errors))
        end
    end

    @testset "Ensembl homology download retries" begin
        response(status::Integer,
            body::AbstractString = "") = Iduna.ThorAxeMSA.HTTP.Response(status, Vector{UInt8}(body))

        attempts = Ref(0)
        data = Iduna.ThorAxeMSA._fetch_ensembl_homology_data(
            "homo_sapiens", "ENSG00000000001.1";
            retries = 3,
            sleep_seconds = 0,
            http_get = (url;
                headers,
                retry,
                status_exception) -> begin
                attempts[] += 1
                @test occursin("/homology/id/homo_sapiens/ENSG00000000001", url)
                @test headers == Iduna.ThorAxeMSA._ENSEMBL_JSON_HEADERS
                @test retry == false
                @test status_exception == false
                attempts[] == 1 ? response(500) : response(200, """{"data": []}""")
            end)
        @test get(data, "data", nothing) !== nothing
        @test attempts[] == 2

        attempts[] = 0
        @test_throws ErrorException Iduna.ThorAxeMSA._fetch_ensembl_homology_data(
            "homo_sapiens", "ENSG00000000001.1";
            retries = 3,
            sleep_seconds = 0,
            http_get = (url; kwargs...) -> begin
                attempts[] += 1
                response(400)
            end)
        @test attempts[] == 1
    end

    @testset "transcript_query receives orthology and effective species list" begin
        mktempdir() do tmp
            captured = Ref{Cmd}()
            runner = command -> (captured[] = command)
            Iduna.ThorAxeMSA._run_transcript_query_once(
                "ENSG00000000001", tmp, "homo_sapiens",
                "homo_sapiens,mus_musculus",
                joinpath(tmp, "stdout.log"), joinpath(tmp, "stderr.log");
                orthology = "1:n",
                runner = runner)
            parts = captured[].exec
            @test "--orthology" in parts
            @test parts[findfirst(==("--orthology"), parts) + 1] == "1:n"
            @test "--specieslist" in parts
            @test parts[findfirst(==("--specieslist"), parts) + 1] ==
                  "homo_sapiens,mus_musculus"

            captured_without_species = Ref{Cmd}()
            Iduna.ThorAxeMSA._run_transcript_query_once(
                "ENSG00000000001", tmp, nothing, nothing,
                joinpath(tmp, "stdout_no_species.log"),
                joinpath(tmp, "stderr_no_species.log");
                orthology = "1:1",
                runner = command -> (captured_without_species[] = command))
            no_species_parts = captured_without_species[].exec
            @test "ENSG00000000001" in no_species_parts
            @test "--orthology" in no_species_parts
        end
    end

    @testset "transcript_query retry helper branches" begin
        mktempdir() do tmp
            gene_core = "ENSG00000000001"
            stdout_log = joinpath(tmp, "stdout.log")
            stderr_log = joinpath(tmp, "stderr.log")
            tmp_gene_dir = joinpath(tmp, gene_core)

            failed_action = Iduna.ThorAxeMSA._transcript_query_attempt_action!(
                gene_core, tmp, "homo_sapiens", "homo_sapiens,mus_musculus",
                stdout_log, stderr_log, tmp_gene_dir;
                attempt = 1,
                attempts = 1,
                orthology = "1:1",
                runner = command -> nothing)
            @test failed_action === :failed

            invalid_action = @test_logs (:warn, r"invalid bundle") Iduna.ThorAxeMSA._transcript_query_attempt_action!(
                gene_core, tmp, "homo_sapiens", "homo_sapiens,mus_musculus",
                stdout_log, stderr_log, tmp_gene_dir;
                attempt = 1,
                attempts = 2,
                orthology = "1:1",
                runner = command -> nothing)
            @test invalid_action === :retry

            retry_calls = Ref(0)
            retried_specieslists = String[]
            slept = Float64[]
            retry_runner = command -> begin
                retry_calls[] += 1
                parts = command.exec
                specieslist_index = findfirst(==("--specieslist"), parts)
                push!(retried_specieslists,
                    specieslist_index === nothing ? "" : parts[specieslist_index + 1])
                retry_calls[] == 2 &&
                    write_test_ensembl_bundle(joinpath(pwd(), gene_core))
                nothing
            end
            @test_logs (:warn, r"invalid bundle") Iduna.ThorAxeMSA._run_transcript_query_with_retries!(
                tmp_gene_dir, gene_core, tmp, "homo_sapiens",
                "homo_sapiens,mus_musculus", 2, stdout_log, stderr_log;
                orthology = "1:1",
                runner = retry_runner,
                sleep_fn = seconds -> push!(slept, seconds))
            @test retry_calls[] == 2
            @test retried_specieslists ==
                  ["homo_sapiens,mus_musculus", "homo_sapiens,mus_musculus"]
            @test slept == [1.0]

            rm(tmp_gene_dir; recursive = true, force = true)
            failed_calls = Ref(0)
            Iduna.ThorAxeMSA._run_transcript_query_with_retries!(
                tmp_gene_dir, gene_core, tmp, "homo_sapiens",
                "homo_sapiens,mus_musculus", 1, stdout_log, stderr_log;
                orthology = "1:1",
                runner = command -> (failed_calls[] += 1),
                sleep_fn = seconds -> error("failed attempts should not sleep"))
            @test failed_calls[] == 1

            target = Iduna.ResolvedTarget(;
                input_id = "ENST00000000001.1",
                input_kind = :ensembl_transcript,
                ensembl_gene_id = "$(gene_core).1",
                transcript_id = "ENST00000000001.1",
                species = "homo_sapiens")
            failed_workdir = joinpath(tmp, "failed_query")
            mkpath(failed_workdir)
            try
                Iduna.ThorAxeMSA._ensure_transcript_query(
                    target, failed_workdir;
                    specieslist = "homo_sapiens,mus_musculus",
                    max_retries = 1,
                    orthology = "1:1",
                    transcript_query_runner = command -> nothing,
                    sleep_fn = seconds -> nothing)
                error("expected transcript_query failure")
            catch err
                @test err isa ErrorException
                @test occursin("transcript_query did not create a valid Ensembl bundle",
                    sprint(showerror, err))
                @test occursin("smaller curated specieslist", sprint(showerror, err))
            end
        end
    end

    @testset "cached ThorAxe MSA result reuse" begin
        mktempdir() do tmp
            source = joinpath(tmp, "cached_input")
            ensembl = joinpath(source, "Ensembl")
            mkpath(ensembl)
            for file in Iduna.ThorAxeMSA._REQUIRED_ENSEMBL_FILES
                write(joinpath(ensembl, file), "$(file)\n")
            end

            workdir = joinpath(tmp, "work")
            target = Iduna.ResolvedTarget(;
                input_id = "ENST00000000001.1",
                input_kind = :ensembl_transcript,
                ensembl_gene_id = "ENSG00000000001.1",
                transcript_id = "ENST00000000001.1")
            input_dir = Iduna.ThorAxeMSA._ensure_transcript_query(target, workdir;
                specieslist = "homo_sapiens",
                cached_input_dir = source,
                orthology = "1:1",
                overwrite = false)
            metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
                input_dir, target, [10.0];
                sample_count = 0,
                sample_fraction = 0.8,
                sample_seed = UInt64(7),
                requested_sample_seed = 7,
                effective_specieslist = "homo_sapiens",
                orthology = "1:1",
                specieslist_filter = true,
                biomart_datasets_filter = true)

            paths = Iduna.ThorAxeMSA._pid_sample_paths(workdir, 10.0, 0)
            write(joinpath(tmp, "seed.fasta"), ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")
            seed_msa = Iduna.ThorAxeMSA.read_file(
                joinpath(tmp, "seed.fasta"), Iduna.ThorAxeMSA.FASTA)
            Iduna.Utils.set_s_exon_annotations!(seed_msa, "00", ['0' => "1_0"])
            mkpath(dirname(paths.fasta_path))
            Iduna.ThorAxeMSA.write_file(paths.fasta_path, seed_msa, Iduna.ThorAxeMSA.FASTA)
            Iduna.ThorAxeMSA.write_file(
                paths.stockholm_path, seed_msa, Iduna.ThorAxeMSA.Stockholm)
            Iduna.ThorAxeMSA._write_candidate_sample_inputs(
                paths, seed_msa, ["homo_sapiens", "mus_musculus"], [1, 2])
            thoraxe_dir = Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(workdir, 10.0, 0)
            mkpath(thoraxe_dir)
            write(joinpath(thoraxe_dir, "path_table.csv"),
                "TranscriptIDCluster,Path\nENST00000000001,start/1_0/stop\n")

            candidate = (;
                msa = seed_msa,
                fasta_path = paths.fasta_path,
                stockholm_path = paths.stockholm_path,
                sequence_fasta = paths.sequence_fasta,
                species_file = paths.species_file,
                workdir)
            row = Iduna.ThorAxeMSA._candidate_summary_row(
                target, candidate, 10.0, 1,
                (; eligible = true, status = "ok", issue = missing);
                sample_count = metadata.pid_sample_count,
                sample_fraction = metadata.pid_sample_fraction,
                sample_seed = metadata.pid_sample_seed,
                metadata)
            summary_path = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir),
                "candidate_summary.csv")
            Iduna.ThorAxeMSA._summarize_candidate_scores([row], summary_path)
            seed = Iduna.SeedSelection(;
                pid = 10.0,
                median_identity = missing,
                mean_identity = missing,
                stockholm_path = paths.stockholm_path,
                fasta_path = paths.fasta_path,
                summary_path)
            Iduna.ThorAxeMSA._mark_selected_candidates!(summary_path, [seed])
            @test Iduna.ThorAxeMSA._cached_selected_seeds(
                summary_path, workdir, metadata) !== nothing

            stale_seed = Iduna.ThorAxeMSA.read_file(
                joinpath(tmp, "seed.fasta"), Iduna.ThorAxeMSA.FASTA)
            Iduna.ThorAxeMSA.write_file(
                paths.stockholm_path, stale_seed, Iduna.ThorAxeMSA.Stockholm)
            @test Iduna.ThorAxeMSA._cached_selected_seeds(
                summary_path, workdir, metadata) === nothing
            Iduna.ThorAxeMSA.write_file(
                paths.stockholm_path, seed_msa, Iduna.ThorAxeMSA.Stockholm)

            result = @test_logs (:info, r"Reusing cached ThorAxe MSA candidates") match_mode=:any Iduna.ThorAxeMSA.build_thoraxe_msa(
                target, workdir;
                pid_thresholds = [10.0],
                specieslist = "homo_sapiens",
                cached_thoraxe_input_dir = source,
                pid_sample_count = 0,
                pid_sample_seed = 7)
            @test result.input_dir == input_dir
            @test result.pid_sample_seed == UInt64(7)
            @test result.pid_sample_count == 0
            @test result.status === :ok
            @test isempty(result.warnings)
            @test only(result.seeds).pid == 10.0
            @test result.baseline_fastas == [paths.fasta_path]
            @test result.baseline_stockholms == [paths.stockholm_path]
            @test result.sequence_fastas == [paths.sequence_fasta]
            @test result.species_files == [paths.species_file]
            @test result.thoraxe_dirs == [thoraxe_dir]

            function write_cached_test_thoraxe_dir(thoraxe_root::AbstractString)
                msa_dir = joinpath(thoraxe_root, "msa")
                phylosofs_dir = joinpath(thoraxe_root, "phylosofs")
                mkpath(msa_dir)
                mkpath(phylosofs_dir)
                write(joinpath(thoraxe_root, "path_table.csv"),
                    "TranscriptIDCluster,Path\nENST00000000001,start/1_0/stop\n")
                write(joinpath(thoraxe_root, "s_exon_table.csv"),
                    "GeneID,Species,TranscriptIDCluster,S_exonID,S_exon_Sequence\n" *
                    "ENSG00000000001,homo_sapiens,ENST00000000001,1_0,AA\n" *
                    "ORTHO1,mus_musculus,ORTHO1,1_0,AB\n")
                write(joinpath(msa_dir, "msa_s_exon_1_0.fasta"),
                    ">ENSG00000000001\nAA\n>ORTHO1\nAB\n")
                write(joinpath(phylosofs_dir, "s_exons.tsv"), "1_0\ta\n")
                write(joinpath(phylosofs_dir, "transcripts.pir"),
                    ">P1;ENSG00000000001 ENST00000000001 a\naa\nAA*\n")
                return thoraxe_root
            end

            scored_workdir = joinpath(tmp, "scored_work")
            scored_input = Iduna.ThorAxeMSA._ensure_transcript_query(
                target, scored_workdir;
                cached_input_dir = source,
                overwrite = false)
            scored_metadata = Iduna.ThorAxeMSA._candidate_run_metadata(
                scored_input, target, [20.0, 30.0];
                sample_count = 0,
                sample_fraction = 0.8,
                sample_seed = UInt64(9),
                requested_sample_seed = 9,
                effective_specieslist = nothing,
                orthology = "1:1",
                specieslist_filter = false,
                biomart_datasets_filter = false)
            scored_summary = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(scored_workdir),
                "candidate_summary.csv")
            mkpath(dirname(scored_summary))
            scored_candidate = (;
                msa = seed_msa,
                fasta_path = "unselected.fasta",
                stockholm_path = "unselected.sto",
                sequence_fasta = "unselected_sequences.fasta",
                species_file = "unselected_species.txt",
                workdir = scored_workdir)
            scored_rows = [
                Iduna.ThorAxeMSA._candidate_summary_row(
                    target, scored_candidate, 20.0, 1,
                    (; eligible = true, status = "ok", issue = missing);
                    sample_count = scored_metadata.pid_sample_count,
                    sample_fraction = scored_metadata.pid_sample_fraction,
                    sample_seed = scored_metadata.pid_sample_seed,
                    metadata = scored_metadata),
                Iduna.ThorAxeMSA._candidate_summary_row(
                    target, scored_candidate, 30.0, 2,
                    (; eligible = true, status = "ok", issue = missing);
                    sample_count = scored_metadata.pid_sample_count,
                    sample_fraction = scored_metadata.pid_sample_fraction,
                    sample_seed = scored_metadata.pid_sample_seed,
                    metadata = scored_metadata)
            ]
            Iduna.ThorAxeMSA._summarize_candidate_scores(scored_rows, scored_summary)
            for pid in (20.0, 30.0)
                write_cached_test_thoraxe_dir(
                    Iduna.ThorAxeMSA._pid_sample_thoraxe_dir(scored_workdir, pid, 0))
            end
            scored_result = Iduna.ThorAxeMSA.build_thoraxe_msa(target, scored_workdir;
                pid_thresholds = [20.0, 30.0],
                cached_thoraxe_input_dir = source,
                specieslist_filter = false,
                biomart_datasets_filter = false,
                pid_sample_count = 0,
                pid_sample_seed = 9)
            @test scored_result.status === :ok
            @test scored_result.pid_sample_count == 0
            @test [seed.pid for seed in scored_result.seeds] == [20.0, 30.0]
            @test length(scored_result.baseline_stockholms) == 2
            @test all(isfile, scored_result.baseline_stockholms)
            @test isfile(scored_result.pid_summary)
        end
    end

    @testset "orphan ThorAxe candidates do not satisfy current stage identity" begin
        mktempdir() do tmp
            workdir = joinpath(tmp, "work")
            paths = Iduna.ThorAxeMSA._pid_sample_paths(workdir, 10.0, 0)
            mkpath(dirname(paths.fasta_path))
            write(paths.fasta_path, ">stale\nAA\n")
            write(paths.stockholm_path, "# STOCKHOLM 1.0\nstale AA\n//\n")
            summary_path = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir),
                "candidate_summary.csv")
            stage_identity = (; target = "current")
            stage_cache = Iduna.ThorAxeMSA._thoraxe_msa_stage_cache(
                workdir, summary_path, (;), stage_identity; overwrite = false)
            @test stage_cache.cache.status === :missing
            @test stage_cache.summary_matches === false
            prepared = Iduna.ThorAxeMSA._prepare_thoraxe_msa_stage!(
                workdir, summary_path, stage_identity, stage_cache; overwrite = false)
            @test prepared.local_artifacts_are_current === false
            @test prepared.force_pid_rerun === true
        end
    end

    @testset "stale ThorAxe manifest cleanup removes old MSA dir" begin
        mktempdir() do tmp
            workdir = joinpath(tmp, "work")
            msa_dir = Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir)
            mkpath(msa_dir)
            write(joinpath(msa_dir, "stale.txt"), "stale")
            summary_path = joinpath(msa_dir, "candidate_summary.csv")
            stage_identity = (; target = "current")
            stage_cache = (;
                cache = (; reusable = true, status = :done, warning = nothing),
                has_manifest = true,
                summary_matches = false)
            prepared = @test_logs (:warn, r"manifest matched") Iduna.ThorAxeMSA._prepare_thoraxe_msa_stage!(
                workdir, summary_path, stage_identity, stage_cache; overwrite = false)
            @test prepared.force_pid_rerun === true
            @test !isfile(joinpath(msa_dir, "stale.txt"))
        end
    end

    @testset "ThorAxe MSA stage failure records manifest" begin
        mktempdir() do tmp
            workdir = joinpath(tmp, "work")
            target = Iduna.ResolvedTarget(;
                input_id = "ENST00000000001.1",
                input_kind = :ensembl_transcript,
                ensembl_gene_id = "ENSG00000000001.1",
                transcript_id = "ENST00000000001.1")
            input_dir = joinpath(tmp, "input")
            write_test_ensembl_bundle(input_dir)
            summary_path = joinpath(Iduna.ThorAxeMSA._thoraxe_msa_dir(workdir),
                "candidate_summary.csv")
            stage_identity = (; target = target.ensembl_gene_id)
            filters = (;
                species_filter = (; warnings = String[]),
                biomart_filter = (; warnings = String[]))
            prepared = (; action = :run, force_pid_rerun = true)
            failing_runner = (args...; kwargs...) -> error("stage boom")
            @test_throws ErrorException Iduna.ThorAxeMSA._run_thoraxe_msa_stage_with_failure_state!(
                failing_runner, target, input_dir, workdir, summary_path, [10.0],
                nothing, (;), stage_identity, filters, prepared;
                pid_sample_count = 0,
                pid_sample_fraction = 0.8,
                sample_seed = UInt64(1))
            state = Iduna.Utils._read_stage_state(
                Iduna.ThorAxeMSA._thoraxe_msa_stage_dir(workdir))
            @test state["status"] == "failed"
            @test state["action"] == "run"
            @test state["exception"]["type"] == "ErrorException"
            @test state["exception"]["message"] == "stage boom"
        end
    end
end
