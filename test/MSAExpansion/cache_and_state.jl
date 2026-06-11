
@test_throws ErrorException Iduna.MSAExpansion.expand_msa(
    target, seed, joinpath(tmp, "invalid_symfrac"); mmseqs_db = db,
    hmmbuild_symfrac = 1.1)
@test !isdir(joinpath(tmp, "invalid_symfrac", "expansion"))

function write_outputs(outputs; centroids::Bool = false)
    mkpath(dirname(outputs.full_stockholm))
    write(outputs.full_stockholm,
        "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
    write(outputs.match_stockholm,
        "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
    write(outputs.a3m_path, "cached a3m\n")
    write(outputs.hits_fasta, ">seed one\nACDE\n>hit one\nACDFG\n>hit_two\nACD\n")
    if centroids
        mkpath(dirname(outputs.centroid_full_stockholm))
        write(outputs.centroid_full_stockholm,
            "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
        write(outputs.centroid_match_stockholm,
            "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed ACDE\n#=GC SExonCode 0000\n//\n")
        write(outputs.centroid_a3m_path, "cached centroid a3m\n")
        write(outputs.centroid_hits_fasta, ">hit one\nACDF\n")
    end
    return outputs
end

function build_self_hit_mmseqs_db(root::AbstractString)
    mmseqs = Iduna.MSAExpansion.MMseqs2_jll.mmseqs()
    mkpath(root)
    fasta = joinpath(root, "self_hit.fasta")
    write(fasta, ">seed\nACDEFGHIKLMNPQRSTVWY\n")
    db = joinpath(root, "self_hit_db")
    seq_db = string(db, "_seq")
    aln_db = string(db, "_aln")
    tmp_dir = joinpath(root, "tmp")
    logs_dir = joinpath(root, "logs")
    mkpath.((tmp_dir, logs_dir))
    Iduna.MSAExpansion._run_labeled(
        `$(mmseqs) createdb $fasta $db`, "createdb", logs_dir)
    Iduna.MSAExpansion._run_labeled(
        `$(mmseqs) createdb $fasta $seq_db`, "createdb_seq", logs_dir)
    Iduna.MSAExpansion._run_labeled(
        `$(mmseqs) search $db $seq_db $aln_db $tmp_dir -a --threads 1`,
        "search", logs_dir)
    return db
end

archive_seed_sto = joinpath(tmp, "archive_seed.sto")
archive_seed_fasta = joinpath(tmp, "archive_seed.fasta")
write(archive_seed_sto,
    "# STOCKHOLM 1.0\n#=GF SExonCodeMap \"0\"=>\"1_0\"\nseed AC\n#=GC SExonCode 00\n//\n")
write(archive_seed_fasta, ">seed\nAC\n")
archive_ctx = (;
    seed_dir = joinpath(tmp, "archive_seeds"),
    seed_stockholm = archive_seed_sto,
    seed_fasta = archive_seed_fasta)
mkpath(archive_ctx.seed_dir)
archived_with_fasta = Iduna.MSAExpansion._archive_expansion_seed(
    seed, archive_ctx)
@test isfile(archived_with_fasta.archived_seed_fasta)
@test read(archived_with_fasta.archived_seed_fasta, String) == ">seed\nAC\n"
@test archived_with_fasta.seed_s_exon_codes == "00"

run_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_10.00")
outputs = write_outputs(
    Iduna.MSAExpansion._expansion_output_paths(run_dir, transcript_id))
identity = Iduna.MSAExpansion._expansion_identity(
    target, seed, seed_sto, nothing, db;
    match_mode = 1,
    match_ratio = nothing,
    hmmbuild_symfrac = 0.0,
    centroids = false)
Iduna.MSAExpansion._write_step_state(run_dir, :done, identity, outputs)
state = Iduna.MSAExpansion._read_step_state(run_dir)
@test state.status == "done"
@test state.step == "msa_expansion"
@test state.stage_key == "expansion:$(gene_id):$(transcript_id):pid_10.00"
@test Iduna.MSAExpansion._step_state_unreadable_message(nothing) ==
      "state file disappeared while reading"
hits_fasta = outputs.hits_fasta
@test_throws MethodError Iduna.MSAExpansion.expand_msa(
    target, seed, tmp; mmseqs_db = db, threads = 1)

cached = @test_logs (:info, r"Reusing cached MSA expansion") match_mode=:any Iduna.MSAExpansion.expand_msa(
    target, seed, tmp; mmseqs_db = db, mmseqs_threads = 1)
cached_with_new_mmseqs_threads = Iduna.MSAExpansion.expand_msa(
    target, seed, tmp; mmseqs_db = db, mmseqs_threads = 8)
hit_counts = Iduna.MSAExpansion._cached_hit_counts(hits_fasta, Set(["seed"]))
@test cached.status === :skipped
@test cached_with_new_mmseqs_threads.status === :skipped
@test cached.n_hits == hit_counts.n_hits
@test cached.n_hits == 3
@test cached.n_new_hits == 2

overwrite_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_15.00")
overwrite_marker = joinpath(overwrite_dir, "old_output.txt")
mkpath(overwrite_dir)
write(overwrite_marker, "stale\n")
overwrite_cache = @test_logs (
    :info, r"Clearing existing MSA expansion run directory") Iduna.MSAExpansion._prepare_expansion_cache!(
    overwrite_dir, tmp, identity, outputs; overwrite = true)
@test overwrite_cache.reusable === false
@test isempty(overwrite_cache.cache_warnings)
@test !ispath(overwrite_dir)

changed_mode_identity = Iduna.MSAExpansion._expansion_identity(
    target, seed, seed_sto, nothing, db;
    match_mode = 3,
    match_ratio = nothing,
    hmmbuild_symfrac = 0.0,
    centroids = false)
changed_mode = Iduna.MSAExpansion._classify_step_state(
    run_dir, changed_mode_identity, outputs)
@test changed_mode.reusable === false
@test changed_mode.status === :stale

changed_ratio_identity = Iduna.MSAExpansion._expansion_identity(
    target, seed, seed_sto, nothing, db;
    match_mode = 1,
    match_ratio = 0.75,
    hmmbuild_symfrac = 0.0,
    centroids = false)
changed_ratio = Iduna.MSAExpansion._classify_step_state(
    run_dir, changed_ratio_identity, outputs)
@test changed_ratio.reusable === false
@test changed_ratio.status === :stale

changed_symfrac_identity = Iduna.MSAExpansion._expansion_identity(
    target, seed, seed_sto, nothing, db;
    match_mode = 1,
    match_ratio = nothing,
    hmmbuild_symfrac = 0.5,
    centroids = false)
changed_symfrac = Iduna.MSAExpansion._classify_step_state(
    run_dir, changed_symfrac_identity, outputs)
@test changed_symfrac.reusable === false
@test changed_symfrac.status === :stale

write(seed_sto, "# STOCKHOLM 1.0\nseed ACDF\n//\n")
changed_seed_identity = Iduna.MSAExpansion._expansion_identity(
    target, seed, seed_sto, nothing, db;
    match_mode = 1,
    match_ratio = nothing,
    hmmbuild_symfrac = 0.0,
    centroids = false)
changed_seed = Iduna.MSAExpansion._classify_step_state(
    run_dir, changed_seed_identity, outputs)
@test changed_seed.reusable === false
@test changed_seed.status === :stale
write(seed_sto, "# STOCKHOLM 1.0\nseed ACDE\n//\n")

centroids_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_20.00")
centroids_seed = Iduna.SeedSelection(;
    pid = 20.0,
    epli = 100.0,
    stockholm_path = seed_sto,
    summary_path = joinpath(tmp, "seed_summary.csv")
)
centroids_outputs = write_outputs(
    Iduna.MSAExpansion._expansion_output_paths(
        centroids_dir, transcript_id; centroids = true);
    centroids = true)
centroids_identity = Iduna.MSAExpansion._expansion_identity(
    target, centroids_seed, seed_sto, nothing, db;
    match_mode = 1,
    match_ratio = nothing,
    hmmbuild_symfrac = 0.0,
    centroids = true)
Iduna.MSAExpansion._write_step_state(
    centroids_dir, :done, centroids_identity, centroids_outputs)
write(hits_fasta, ">seed one\nACDE\n>hit one\nACDF\n>hit_two\nACDG\n")

cached_centroids = Iduna.MSAExpansion.expand_msa(
    target, centroids_seed, tmp; mmseqs_db = db, centroids = true)
@test cached_centroids.status === :skipped
@test isfile(centroids_outputs.s_exon_blocks_tsv)
@test isfile(centroids_outputs.centroid_s_exon_blocks_tsv)

missing_centroid_dir = joinpath(
    tmp, "expansion", gene_id, transcript_id, "pid_30.00")
missing_centroid_outputs = write_outputs(
    Iduna.MSAExpansion._expansion_output_paths(
    missing_centroid_dir, transcript_id; centroids = true))
missing_centroid_seed = Iduna.SeedSelection(;
    pid = 30.0,
    epli = 100.0,
    stockholm_path = seed_sto,
    summary_path = joinpath(tmp, "seed_summary.csv")
)
missing_centroid_identity = Iduna.MSAExpansion._expansion_identity(
    target, missing_centroid_seed, seed_sto, nothing, db;
    match_mode = 1,
    match_ratio = nothing,
    hmmbuild_symfrac = 0.0,
    centroids = true)
Iduna.MSAExpansion._write_step_state(
    missing_centroid_dir, :done, missing_centroid_identity,
    missing_centroid_outputs)
missing_centroids = Iduna.MSAExpansion._classify_step_state(
    missing_centroid_dir, missing_centroid_identity, missing_centroid_outputs)
@test missing_centroids.reusable === false
@test missing_centroids.status === :unfinished

incomplete_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_35.00")
incomplete_outputs = Iduna.MSAExpansion._expansion_output_paths(
    incomplete_dir, transcript_id)
mkpath(incomplete_dir)
incomplete = Iduna.MSAExpansion._classify_step_state(
    incomplete_dir, identity, incomplete_outputs)
@test incomplete.reusable === false
@test incomplete.status === :unfinished
@test occursin("incomplete outputs", incomplete.warning)

missing_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_36.00")
missing_outputs = Iduna.MSAExpansion._expansion_output_paths(
    missing_dir, transcript_id)
missing_state = Iduna.MSAExpansion._classify_step_state(
    missing_dir, identity, missing_outputs)
@test missing_state.reusable === false
@test missing_state.status === :missing
@test missing_state.warning === nothing

unreadable_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_37.00")
unreadable_outputs = Iduna.MSAExpansion._expansion_output_paths(
    unreadable_dir, transcript_id)
mkpath(unreadable_dir)
write(Iduna.MSAExpansion._step_state_path(unreadable_dir), "{not json")
unreadable = Iduna.MSAExpansion._classify_step_state(
    unreadable_dir, identity, unreadable_outputs)
@test unreadable.reusable === false
@test unreadable.status === :stale
@test occursin("Could not read MSA expansion stage_state.json", unreadable.warning)

legacy_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_40.00")
legacy_outputs = write_outputs(
    Iduna.MSAExpansion._expansion_output_paths(legacy_dir, transcript_id))
legacy_seed = Iduna.SeedSelection(;
    pid = 40.0,
    epli = 100.0,
    stockholm_path = seed_sto,
    summary_path = joinpath(tmp, "seed_summary.csv")
)
legacy_identity = Iduna.MSAExpansion._expansion_identity(
    target, legacy_seed, seed_sto, nothing, db;
    match_mode = 1,
    match_ratio = nothing,
    hmmbuild_symfrac = 0.0,
    centroids = false)
legacy = Iduna.MSAExpansion._classify_step_state(
    legacy_dir, legacy_identity, legacy_outputs)
@test legacy.reusable === false
@test legacy.status === :stale
@test occursin("no stage_state.json", legacy.warning)

stale_failure_dir = joinpath(tmp, "expansion", gene_id, transcript_id, "pid_50.00")
stale_failure_outputs = write_outputs(
    Iduna.MSAExpansion._expansion_output_paths(stale_failure_dir, transcript_id))
