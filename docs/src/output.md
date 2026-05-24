# Output Layout

Iduna keeps ThorAxe's own output layout untouched and writes package metadata
around it. By default, `workdir` is a directory named after the input ID.

```text
<workdir>/
  target.json
  result.json
  logs/
    thoraxe/
  sequences/
    uniprot/
    ensembl_proteins/
  thoraxe_input/
    Ensembl/
  thoraxe_msa/
    runs/
      pid_10.00/
        full/
          thoraxe/
    candidates/
      pid_10.00/
        candidate_msa_full.fasta
        candidate_msa_full.sto
        scores.csv
        sequences/
          candidate_sequences_full.fasta
          candidate_sequences_species_subset_001.fasta
        species/
          candidate_species_full.txt
          candidate_species_subset_001.txt
    candidate_summary.csv
  expansion/  # absent when no_expansion=true / --no-expansion
    <gene>/<transcript>/
      pid_10.00/
        step_state.json
        seeds/
        dbs/
        hmm/
        logs/
        expanded_msa/
        centroid_msa/  # optional, written only with centroids=true / --centroids
  validation/
    pid_10.00/
      stats.csv
      query_vs_uniprot_alignment.txt
```

The output layout is part of the API. `overwrite=true` rebuilds package-owned
subdirectories, but Iduna does not delete files outside the selected work
directory.

After the work directory is created, Iduna always attempts to write
`result.json`. Successful runs keep the returned [`IdunaResult`](@ref) summary
there. Failed runs still throw the original exception, but leave a
machine-readable `result.json` with `status: "error"`, `failed_stage`, any
available partial `target` metadata, accumulated warnings, and exception
details. Timeout failures from logged ThorAxe commands also include the command
and stdout/stderr log paths.

When `thoraxe_input_dir` is supplied, that source bundle is treated as read-only
and copied into `thoraxe_input/`; the copied layout is then preserved like a
fresh `transcript_query` result.

`thoraxe_msa/runs/` keeps the raw ThorAxe output for each PID candidate and its
PID-specific species samples. `thoraxe_msa/candidates/` stores each PID's
full-species candidate `msa_0`, sampled species lists, gap-free sequence files,
and per-PID score CSV. `thoraxe_msa/candidate_summary.csv` records validation
status, eligibility, sample identity statistics, candidate size, and the selected
seed rows.

MSA paths are also available from the returned [`IdunaResult`](@ref):

```julia
result.thoraxe_msa.baseline_stockholms[1]
result.thoraxe_msa.baseline_fastas[1]
result.thoraxe_msa.seeds[1].stockholm_path
result.thoraxe_msa.seeds[1].fasta_path
result.thoraxe_msa.pid_sample_count
result.thoraxe_msa.pid_sample_fraction
result.thoraxe_msa.pid_sample_seed

# Full expansion runs only:
result.expansions[1].match_stockholm
result.expansions[1].a3m_path
```

`result.workdir` is stored as an absolute path. Artifact paths inside
`workdir` are stored relative to it, so use `joinpath(result.workdir, path)` to
open them directly. Paths outside `workdir` remain unchanged.

When `no_expansion=true` or `--no-expansion` is used, Iduna stops after the
ThorAxe MSA stage. `isempty(result.expansions)`, `result.json` contains an empty
`"expansions"` list, and no `expansion/` directory is written. Validation still
writes seed statistics to `validation/pid_<value>/stats.csv`, with expanded-MSA
fields left missing.

By default, Iduna selects one PID seed. Set `pid_sample_count=0` to skip seed
selection and carry every eligible PID candidate forward, producing one
validation directory and, unless `no_expansion=true`, one expansion directory
per selected PID.

`expanded_msa/` is the full MMseqs2 cluster-expanded MSA and remains the main
Iduna result. When `centroids=true` or `--centroids` is used, Iduna also writes
`centroid_msa/` with the MSA built from centroid/consensus hits before cluster
expansion. This side output reuses the same MMseqs2 search result, does not run
a second search, and does not change validation or `result.json`. Its files are
`<transcript>_centroids_full.sto`, `<transcript>_centroids_matchonly.sto`,
`<transcript>_centroids.a3m`, and `<transcript>_centroid_hits_raw.fasta`.

Each expansion PID directory also has `step_state.json`, which records the
input identity used for safe cache reuse. If the seed MSA, expansion parameters,
or centroid request change, Iduna treats the cached expansion as outdated,
warns, and rebuilds that PID directory.
