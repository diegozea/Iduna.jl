## Iduna.jl Release Notes

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
