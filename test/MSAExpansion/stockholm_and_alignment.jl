source = joinpath(tmp, "in.sto")
dest = joinpath(tmp, "out.sto")
write(source, "seq1 AC\n")
Iduna.MSAExpansion.prepare_stockholm_for_mmseqs(source, dest)
text = read(dest, String)
@test startswith(text, "# STOCKHOLM 1.0")
@test endswith(chomp(text), "//")

complete = joinpath(tmp, "complete.sto")
complete_out = joinpath(tmp, "complete_out.sto")
write(complete, "# STOCKHOLM 1.0\nseq1 AC\n//\n")
Iduna.MSAExpansion.prepare_stockholm_for_mmseqs(complete, complete_out)
complete_text = read(complete_out, String)
@test complete_text == "# STOCKHOLM 1.0\nseq1 AC\n//\n"

annotated_seed = joinpath(tmp, "annotated_seed.sto")
sanitized_seed = joinpath(tmp, "sanitized_seed.sto")
write(annotated_seed, """
# STOCKHOLM 1.0
#=GF SExonCodeMap "a"=>"1_0"
seq1 AC
#=GC SExonCode aa
//
""")
Iduna.MSAExpansion.prepare_stockholm_for_mmseqs(annotated_seed, sanitized_seed)
sanitized_text = read(sanitized_seed, String)
@test !occursin("SExonCode", sanitized_text)
@test occursin("seq1 AC", sanitized_text)

empty_sto = joinpath(tmp, "empty.sto")
touch(empty_sto)
@test Iduna.MSAExpansion.normalize_stockholm_annotations!(empty_sto) == empty_sto
@test isempty(read(empty_sto, String))

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

noisy_hits_tsv = joinpath(tmp, "noisy_hits.tsv")
write(noisy_hits_tsv, """
too-short
query\t\tACDE
query\thit one\t---
query\thit one\tACD-
query\thit one\tACDE
query\tseed one\tACDF
query\tnew_hit\tacdg
""")
noisy_all_hits,
noisy_filtered_hits = Iduna.MSAExpansion._collect_hits(noisy_hits_tsv, Set(["seed"]))
@test noisy_all_hits ==
      [("hit", "ACD"), ("seed", "ACDF"), ("new_hit", "ACDG")]
@test noisy_filtered_hits == [("hit", "ACD"), ("new_hit", "ACDG")]

logs_dir = joinpath(tmp, "run_logs")
run_labeled_logs,
_ = Test.collect_test_logs() do
    @test Iduna.MSAExpansion._run_labeled(
        `sh -c "printf msa-out; printf msa-err >&2"`, "mock",
        logs_dir) === nothing
end
@test isempty([log
               for log in run_labeled_logs
               if log.message == "Running MSA expansion command."])
@test read(joinpath(logs_dir, "mock_stdout.log"), String) == "msa-out"
@test read(joinpath(logs_dir, "mock_stderr.log"), String) == "msa-err"

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

sexon_sto = joinpath(tmp, "sexon_project.sto")
write(sexon_sto, """
# STOCKHOLM 1.0
seed ACdeFG
hit AC-eFG
#=GC RF xxxxxx
//
""")
sexon_msa = Iduna.MSAExpansion.read_file(
    sexon_sto, Iduna.MSAExpansion.Stockholm; keepinserts = true)
@test Iduna.MSAExpansion.getannotcolumn(sexon_msa, "Aligned", "") == "110011"
archived = (;
    seed_s_exon_codes = "0β23",
    seed_match_s_exon_codes = "0β23",
    seed_s_exon_code_map = [
        '0' => "1_0", 'β' => "2_0", '2' => "3_0", '3' => "4_0"]
)
Iduna.MSAExpansion._restore_s_exon_annotations!(sexon_msa, archived)
@test Iduna.Utils.s_exon_codes(sexon_msa) == "0β..23"
@test Iduna.Utils.s_exon_code_map(sexon_msa)['3'] == "4_0"

annotated_seed_sto = joinpath(tmp, "annotated_seed_rf.sto")
write(annotated_seed_sto, """
# STOCKHOLM 1.0
seed ACdE
#=GC RF xxxx
//
""")
annotated_seed = Iduna.MSAExpansion.read_file(
    annotated_seed_sto, Iduna.MSAExpansion.Stockholm; keepinserts = true)
@test Iduna.MSAExpansion.getannotcolumn(annotated_seed, "Aligned", "") == "1101"
@test Iduna.MSAExpansion._seed_match_s_exon_codes("0123", annotated_seed) ==
      "013"
masked_archived = (;
    seed_s_exon_codes = "0123",
    seed_match_s_exon_codes = "013",
    seed_s_exon_code_map = ['0' => "1_0", '1' => "2_0", '3' => "4_0"]
)
Iduna.MSAExpansion._restore_s_exon_annotations!(sexon_msa, masked_archived)
@test Iduna.Utils.s_exon_codes(sexon_msa) == "01..3."

@test Iduna.MSAExpansion._rf_match_state_mask("xx..xx", 6) ==
      [true, true, false, false, true, true]
@test Iduna.MSAExpansion._rf_match_state_mask("", 6) === nothing
@test Iduna.MSAExpansion._aligned_match_state_mask("", 6) === nothing
plain_sto = joinpath(tmp, "plain_no_annotations.sto")
write(plain_sto, "# STOCKHOLM 1.0\nseed AC\n//\n")
plain_msa = Iduna.MSAExpansion.read_file(
    plain_sto, Iduna.MSAExpansion.Stockholm; keepinserts = true)
@test Iduna.MSAExpansion._match_state_mask(
    plain_msa; default_aligned = true) == [true, true]
plain_msa_without_insert_annotation = Iduna.MSAExpansion.read_file(
    plain_sto, Iduna.MSAExpansion.Stockholm)
@test Iduna.MSAExpansion._match_state_mask(
    plain_msa_without_insert_annotation; default_aligned = false) ==
      [false, false]

gene_id = "ENSG00000198821"
transcript_id = "ENST00000362089.10"
seed_sto = joinpath(tmp, "seed.sto")
write(seed_sto, "# STOCKHOLM 1.0\nseed ACDE\n//\n")
seed = Iduna.SeedSelection(;
    pid = 10.0,
    epli = 100.0,
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
