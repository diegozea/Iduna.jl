## Iduna.jl Release Notes

### Changes from v0.7.0 to v0.8.0

Breaking changes:

- Changed the default ThorAxe `specieslist` to `"ases"`, matching the Ases
  webserver 12-species list. Use `specieslist="all"` or `specieslist=""` to
  restore unrestricted ThorAxe species selection.
- Changed the default PID seed sampling strategy from PID-local sampling to
  `sampling_strategy=:common`, which draws one shared set of species samples
  from the species common to all eligible PID candidates. Use
  `sampling_strategy=:independent` to restore the v0.7 behavior, or
  `sampling_strategy=:input` to sample from the effective input species list.
- Replaced expansion `step_state.json` files with the shared
  `stage_state.json` manifest format. Code that reads expansion state files
  directly should use the new filename and fields.
- Saved result metadata is now move-safe and omits the absolute `workdir`; paths
  for artifacts under the work directory are stored relative to it. The in-memory
  `IdunaResult.workdir` remains absolute.
- `IdunaResult.expansions` can contain `nothing` for a seed whose expansion is
  unavailable, while still remaining indexed like `result.thoraxe_msa.seeds`.

Added:

- Added manifest-backed resumability for pipeline stages. `target`,
  `thoraxe_input`, `thoraxe_msa`, expansion, validation, and `result` stages now
  record `stage_state.json` manifests with input identities, outputs, actions,
  warnings, and failures.
- `overwrite=false` reruns now reuse completed matching stages and rebuild
  missing, stale, failed, or incomplete stages. Copied or generated
  `transcript_query` bundles are fingerprinted and reused when their target,
  species list, orthology, and filter options still match.
- Added the public `load_result(workdir)` API to reconstruct an `IdunaResult`
  from current result artifacts without rerunning pipeline stages or writing
  files, including after the result directory has been moved or copied.
- Added `sampling_strategy` to the Julia API and `--sampling-strategy` to the
  app. Supported strategies are `common`, `independent`, and `input`.
- Added shared PID sample species files under `thoraxe_msa/samples/species/`
  for shared sampling strategies, plus `sampling_strategy` metadata in
  ThorAxe MSA results and candidate summaries.
- Parallelized ThorAxe PID sample scoring across Julia threads. Start Julia
  with `--threads` or set `JULIA_NUM_THREADS` to speed up the repeated ThorAxe
  sample runs; `--mmseqs-threads` still controls only MMseqs2 expansion.
- Improved `IdunaResult` text display with a compact summary of pipeline
  status, selected seeds, expansion slots, validations, warnings, and errors.
- Expanded progress logging for ThorAxe transcript_query preparation, PID
  candidate generation, sample scoring, cache reuse, and retry/failure paths.

Internal changes:

- Kept HTTP.jl v1 compatibility because HTTP.jl v2 interrupts
  JuliaRegistrator AutoMerge, and added explicit CodecZlib-based gzip response
  decoding.
- Added tests for stage resumability, move-safe result loading, partial result
  round-tripping, sampling strategies, threaded PID scoring, display output, and
  fallback/error paths.
- Updated CI and documentation dependencies, including `codecov/codecov-action`
  v6, `julia-actions/cache` v3, and Documenter compatibility for 1.17.0.
- Added the WIP repository status badge.

### Changes from v0.6.0 to v0.7.0

Breaking changes:

- Removed the transcript_query and ThorAxe wall-clock limit controls from the
  API and CLI: `transcript_query_timeout_seconds`,
  `transcript_query_timeout_max_seconds`, `allow_specieslist_timeout_fallback`,
  and `thoraxe_timeout_seconds`. Use an external scheduler or process wrapper
  when a hard runtime limit is needed.
- transcript_query retries now preserve the same effective species list. If an
  invalid or incomplete bundle persists, Iduna fails explicitly and suggests
  trying a smaller curated `specieslist`.

Added:

- Added `@info` progress logging for the main Iduna pipeline stages, including
  output preparation, target resolution, ThorAxe MSA generation, MSA expansion,
  validation, and result writing.
- Added lower-level progress logging for ThorAxe species filtering, cache use,
  PID scoring, seed selection, MMseqs and HMMER work, and expansion outputs.

Internal changes:

- Updated tests for the timeout removal and progress logging behavior.
- Updated documentation and examples to remove references to the removed
  timeout controls.

### Changes from v0.5.0 to v0.6.0

Breaking changes:

- Iduna now needs ThorAxe runs to contain the extra files used to track s-exon
  columns. Older cached ThorAxe outputs that do not have those files are rebuilt
  or completed instead of being reused as they are.
- `result.json` and result summaries now include paths to the new s-exon block
  tables. Programs that expect an exact list of JSON fields should allow these
  new fields.
- Seed and expansion result objects now include `s_exon_blocks_tsv`, the path to
  the s-exon block table. Code that checks the exact fields of these objects may
  need to be updated.

Added:

- Added a way to keep track of which ThorAxe s-exon each alignment column came
  from. This information is stored in Stockholm files as MIToS annotations:
  `SExonCode` has one short symbol per column, and `SExonCodeMap` explains
  which s-exon ID each symbol means.
- Added `*_s_exon_blocks.tsv` files for seed, expanded, and centroid
  alignments. These tables group neighboring columns that come from the same
  s-exon, which makes them easier to use for plots and manual inspection.
- Iduna keeps the s-exon column labels when the seed MSA is expanded with HMMER
  and MMseqs. Insert columns that cannot be assigned to one s-exon are left
  unassigned.
- Iduna now reads ThorAxe's PhyloSofS s-exon symbols and `s_exons.tsv` mapping,
  including `0_` s-exons when they have protein sequence.
- Added documentation showing how to read the block tables and the MIToS
  annotations.

Internal changes:

- Replaced JSON3 and CodecZlib usage with JSON.jl and HTTP body decoding.
- Added tests for ThorAxe retry, cache, scoring, and s-exon annotation paths
  that do not depend on live web services.
- Added CodeComplexity checks to the test suite.
- Skipped the live MAPK8 integration smoke test on GitHub Actions while keeping
  it enabled for local test runs.

### Changes from v0.4.0 to v0.5.0

Breaking changes:

- Changed expansion cache reuse to require a matching `step_state.json` input
  identity. Expansion directories created by older Iduna versions do not have
  this file, so they are treated as outdated and rebuilt once instead of being
  reused solely because output files exist.

Added:

- Added `step_state.json` to each expansion PID directory, recording the seed
  MSA identity, expansion parameters, MMseqs2 database path, requested centroid
  output state, expected output paths, status, warnings, and failed-run
  exception summaries.
- Added checks that rebuild cached expansion outputs when the seed MSA,
  expansion options, MMseqs2 database path, or centroid request differs from the
  recorded expansion state.

### Changes from v0.3.0 to v0.4.0

Breaking changes:

- Replaced singular top-level result fields `IdunaResult.expansion` and
  `IdunaResult.validation` with `IdunaResult.expansions` and
  `IdunaResult.validations`, allowing one result entry per selected PID seed.
- Replaced singular ThorAxe MSA fields such as
  `ThorAxeMSAResult.best_seed`, `baseline_stockholm`, `baseline_fasta`,
  `sequence_fasta`, `species_file`, and `thoraxe_dir` with vector-based fields
  such as `seeds`, `baseline_stockholms`, `baseline_fastas`,
  `sequence_fastas`, `species_files`, and `thoraxe_dirs`.
- Changed `no_expansion=true` results to use empty `expansions` instead of
  `expansion === nothing`; JSON summaries now store `"expansions": []`.
- Moved expansion and validation artifacts into PID-specific directories.

Added:

- Added PID candidate sampling and seed selection controlled by
  `pid_sample_count`, `pid_sample_fraction`, and `pid_sample_seed`.
- Added reproducible PID sampling when `pid_sample_seed` is supplied; otherwise
  a random seed is recorded in the result metadata.
- Added candidate validation that excludes PID candidates with indels relative
  to UniProt and reports substitution-only differences as warnings.
- Added `pid_sample_count=0` mode to skip seed selection and carry every
  eligible PID candidate forward through expansion and validation.
- Added `pid` and `index` selectors to `load_seed_msa` and
  `load_expanded_msa` for results with multiple selected seeds.
- Added `thoraxe_msa/candidate_summary.csv` with candidate eligibility,
  identity scores, sampling metadata, and selected rows.
- Added transcript-query metadata and bundle fingerprint checks for safer
  reuse of cached ThorAxe input.

Internal changes:

- Added `Random`, `SHA`, and `StatsBase` dependencies for PID sampling,
  fingerprinting, and species subset sampling.

### Changes from v0.2.0 to v0.3.0

Breaking changes:

- Changed returned result artifact paths under `workdir` to be relative while
  keeping `IdunaResult.workdir` absolute, making printed result objects and
  JSON summaries more compact. We have kept lower-level readers and validators 
  compatible with returned relative paths by resolving them against their 
  associated work directory.

### Changes from v0.1.0 to v0.2.0

- Added `no_expansion=true` and `--no-expansion` to stop after the ThorAxe MSA
  stage without requiring an MMseqs2 database; in that mode,
  `IdunaResult.expansion` is `nothing` and the JSON summary stores
  `"expansion": null`.
- Added `centroids=true` and `--centroids` to write the centroid-level MSA side
  output before MMseqs2 expands hits to all cluster members.
- Added one-field-per-line pretty printing for `IdunaResult` and the nested
  result objects returned by the pipeline.
