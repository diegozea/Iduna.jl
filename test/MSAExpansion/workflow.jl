stale_failure_seed = Iduna.SeedSelection(;
    pid = 50.0,
    epli = 100.0,
    stockholm_path = seed_sto,
    summary_path = joinpath(tmp, "seed_summary.csv")
)
@test_throws Base.ProcessFailedException Iduna.MSAExpansion.expand_msa(
    target, stale_failure_seed, tmp; mmseqs_db = db, mmseqs_threads = 1)
stale_failure_state = Iduna.MSAExpansion._read_step_state(stale_failure_dir)
@test stale_failure_state.status == "failed"
@test occursin("no stage_state.json", stale_failure_state.warnings[1])
@test occursin("ProcessFailedException", stale_failure_state.exception.type)
cache_warning = read(joinpath(stale_failure_dir, "logs", "cache_warning.log"),
    String)
@test occursin("no stage_state.json", cache_warning)
@test !isfile(stale_failure_outputs.full_stockholm)

success_seed_sto = joinpath(tmp, "success_seed.sto")
write(success_seed_sto, "# STOCKHOLM 1.0\nseed ACDEFGHIKLMNPQRSTVWY\n//\n")
success_seed = Iduna.SeedSelection(;
    pid = 60.0,
    epli = 100.0,
    stockholm_path = success_seed_sto,
    summary_path = joinpath(tmp, "seed_summary.csv")
)
success_db = build_self_hit_mmseqs_db(joinpath(tmp, "self_hit_mmseqs"))
success_logs,
success = Test.collect_test_logs() do
    Iduna.MSAExpansion.expand_msa(
        target, success_seed, tmp; mmseqs_db = success_db, centroids = true,
        mmseqs_threads = 1)
end
@test isempty([log
               for log in success_logs
               if log.message in ("Preparing MSA expansion directories.",
    "Archiving MSA expansion seed.",
    "Converting MMseqs hits.",
    "Writing MSA expansion hits.",
    "Building seed HMM.",
    "Aligning MSA expansion hits.",
    "Writing expanded alignment outputs.",
    "Writing centroid MSA.")])
preparing_log = only([log
                      for log in success_logs
                      if log.message == "Preparing MSA expansion."])
preparing_kwargs = Dict(preparing_log.kwargs)
@test preparing_kwargs[:pid] == 60.0
@test haskey(preparing_kwargs, :mmseqs_db)
@test haskey(preparing_kwargs, :centroids)
@test !haskey(preparing_kwargs, :mmseqs_threads)
workflow_log = only([log
                     for log in success_logs
                     if log.message == "Running MMseqs expansion workflow."])
workflow_kwargs = Dict(workflow_log.kwargs)
@test workflow_kwargs[:mmseqs_threads] == 1
@test !haskey(workflow_kwargs, :pid)
@test !haskey(workflow_kwargs, :tmp_dir)
completion_log = only([log
                       for log in success_logs
                       if log.message == "MSA expansion completed."])
completion_kwargs = Dict(completion_log.kwargs)
@test completion_kwargs[:pid] == 60.0
@test completion_kwargs[:n_new_hits] == 0
@test !haskey(completion_kwargs, :run_dir)
@test success.status === :ok
@test success.n_new_hits == 0
@test isfile(joinpath(success.run_dir, "centroid_msa",
    "$(transcript_id)_centroids.a3m"))
success_state = Iduna.MSAExpansion._read_step_state(success.run_dir)
@test success_state.status == "done"

Iduna.MSAExpansion._write_step_state(legacy_dir, :failed, legacy_identity,
    legacy_outputs; exception = (; type = "ErrorException", message = "failed"))
failed = Iduna.MSAExpansion._classify_step_state(
    legacy_dir, legacy_identity, legacy_outputs)
@test failed.reusable === false
@test failed.status === :failed

Iduna.MSAExpansion._write_step_state(legacy_dir, :unfinished, legacy_identity,
    legacy_outputs)
unfinished = Iduna.MSAExpansion._classify_step_state(
    legacy_dir, legacy_identity, legacy_outputs)
@test unfinished.reusable === false
@test unfinished.status === :unfinished
