## Iduna.jl Release Notes

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
- Removed the transcript_query and ThorAxe wall-clock limit controls from the
  API and CLI: `transcript_query_timeout_seconds`,
  `transcript_query_timeout_max_seconds`, `allow_specieslist_timeout_fallback`,
  and `thoraxe_timeout_seconds`. Use an external scheduler or process wrapper
  when a hard runtime limit is needed.
- transcript_query retries now preserve the same effective species list. If an
  invalid or incomplete bundle persists, Iduna fails explicitly and suggests
  trying a smaller curated `specieslist`.

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
