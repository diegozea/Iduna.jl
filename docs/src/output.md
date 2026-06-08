# Output Layout

Iduna keeps ThorAxe's own output layout untouched and writes package metadata
around it. By default, `workdir` is a directory named after the input ID.

```text
<workdir>/
  target.json
  result.json
  .iduna/
    stages/
      target/
        stage_state.json
      thoraxe_input/
        stage_state.json
      thoraxe_msa/
        stage_state.json
      result/
        stage_state.json
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
    samples/
      species/
        candidate_species_subset_001.txt
    candidates/
      pid_10.00/
        candidate_msa_full.fasta
        candidate_msa_full.sto
        scores.csv
        sequences/
          candidate_sequences_full.fasta
          candidate_sequences_species_subset_001.fasta  # independent sampling only
        species/
          candidate_species_full.txt
          candidate_species_subset_001.txt  # symlink for shared sampling
    candidate_summary.csv
  expansion/  # absent when no_expansion=true / --no-expansion
    <gene>/<transcript>/
      pid_10.00/
        stage_state.json
        seeds/
        dbs/
        hmm/
        logs/
        expanded_msa/
        centroid_msa/  # optional, written only with centroids=true / --centroids
  validation/
    pid_10.00/
      stage_state.json
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
details.

When `thoraxe_input_dir` is supplied, that source bundle is treated as read-only
and copied into `thoraxe_input/`; the copied layout is then preserved like a
fresh `transcript_query` result. Iduna fingerprints the required Ensembl bundle
files and records a `thoraxe_input` stage manifest, so later reruns reuse the
copied or previously generated `transcript_query` result when the target,
species list, orthology, and filter options still match.

Each `stage_state.json` records the stage name, stage key, status, action,
identity hash, declared outputs, warnings, and exception summary when a stage
fails. `result.json` includes a compact `"stages"` list copied from these
manifests, making reruns auditable without using it as the source of truth.

`thoraxe_msa/runs/` keeps the raw ThorAxe output for each PID candidate and its
species samples. `thoraxe_msa/candidates/` stores each PID's full-species
candidate `msa_0`, sampled species list paths, gap-free sequence files, and
per-PID score CSV. With `sampling_strategy=:common` or `:input`, shared sampled
species lists live in `thoraxe_msa/samples/species/` and each PID-local sampled
species path is a symlink to the matching shared file.
`thoraxe_msa/candidate_summary.csv` records validation status, eligibility,
sample identity statistics, candidate size, and the selected seed rows.

MSA paths are also available from the returned [`IdunaResult`](@ref):

```julia
result.thoraxe_msa.baseline_stockholms[1]
result.thoraxe_msa.baseline_fastas[1]
result.thoraxe_msa.seeds[1].stockholm_path
result.thoraxe_msa.seeds[1].fasta_path
result.thoraxe_msa.seeds[1].s_exon_blocks_tsv
result.thoraxe_msa.pid_sample_count
result.thoraxe_msa.pid_sample_fraction
result.thoraxe_msa.pid_sample_seed
result.thoraxe_msa.sampling_strategy

# Full expansion runs only:
expansion = result.expansions[1]
ismissing(expansion) || expansion.match_stockholm
ismissing(expansion) || expansion.s_exon_blocks_tsv
ismissing(expansion) || expansion.a3m_path
```

`result.workdir` is absolute in memory. Saved metadata, including
`result.json`, `target.json`, `candidate_summary.csv`, and stage state files,
stores paths inside `workdir` relative to that directory so the result folder can
be moved or copied. Paths outside `workdir` remain unchanged.

When `no_expansion=true` or `--no-expansion` is used, Iduna stops after the
ThorAxe MSA stage. `isempty(result.expansions)`, `result.json` contains an empty
`"expansions"` list, and no `expansion/` directory is written. Validation still
writes seed statistics to `validation/pid_<value>/stats.csv`, with expanded-MSA
fields left missing.

By default, Iduna selects one PID seed. Set `pid_sample_count=0` to skip seed
selection and carry every eligible PID candidate forward, producing one
validation directory and, unless `no_expansion=true`, one expansion directory
per selected PID. `result.expansions` is seed-indexed and may contain `nothing`
for a seed whose expansion is unavailable.

`expanded_msa/` is the full MMseqs2 cluster-expanded MSA and remains the main
Iduna result. When `centroids=true` or `--centroids` is used, Iduna also writes
`centroid_msa/` with the MSA built from centroid/consensus hits before cluster
expansion. This side output reuses the same MMseqs2 search result, does not run
a second search, and does not change validation or `result.json`. Its files are
`<transcript>_centroids_full.sto`, `<transcript>_centroids_matchonly.sto`,
`<transcript>_centroids.a3m`, `<transcript>_centroids_s_exon_blocks.tsv`, and
`<transcript>_centroid_hits_raw.fasta`.

## S-Exon Column Annotations

Iduna keeps track of which original ThorAxe s-exon each MSA column came from.
This is useful when you want to color the alignment by s-exon blocks, compare
blocks between the seed and expanded MSA, or inspect which part of the original
ThorAxe path supports a region of interest.

For most users, the easiest file to use is the block table:

```julia
seed_blocks = joinpath(
    result.workdir,
    result.thoraxe_msa.seeds[1].s_exon_blocks_tsv,
)

expanded_blocks = joinpath(
    result.workdir,
    something(result.expansions[1]).s_exon_blocks_tsv,
)
```

Each `*_s_exon_blocks.tsv` file has one row per continuous block of columns
with the same s-exon ID:

| column | meaning |
|:---|:---|
| `alignment` | which alignment the row describes, for example `seed`, `match`, or `full` |
| `pid` | the ThorAxe PID threshold used for that seed |
| `code` | the short one-character symbol stored in the MSA |
| `s_exon_id` | the original ThorAxe s-exon ID, such as `12_2` or `0_1` |
| `start_col` | first MSA column in the block, using 1-based numbering |
| `end_col` | last MSA column in the block, using 1-based numbering |
| `n_columns` | number of MSA columns in the block |

The Stockholm MSA files also keep the same information as MIToS annotations.
`SExonCode` is one short symbol per MSA column. `SExonCodeMap` is the key that
translates each symbol back to a ThorAxe s-exon ID.

```text
#=GF SExonCodeMap "a"=>"1_0","b"=>"12_2"
#=GC SExonCode    aaabbb
```

In this small example, columns 1 to 3 come from s-exon `1_0`, and columns 4 to
6 come from s-exon `12_2`.

If you work directly with MIToS, load Stockholm files with `keepinserts=true` so
insert columns are kept:

```julia
using MIToS.MSA

msa = read_file(stockholm_path, Stockholm; keepinserts=true)
column_codes = getannotcolumn(msa, "SExonCode")
code_key = getannotfile(msa, "SExonCodeMap")
```

The TSV block table is usually simpler for plotting, because it already groups
neighboring columns that belong to the same s-exon.

Some details are important when reading these annotations:

- `.` in `SExonCode` means Iduna does not assign that column to an s-exon. This
  usually marks an inserted or unmatched region.
- `-` gap characters inside an s-exon block still keep that s-exon's annotation,
  because the column belongs to that s-exon alignment.
- ThorAxe IDs that start with `0_` are real ThorAxe s-exon IDs. They often mean
  that ThorAxe found a segment that did not align with the other s-exons in the
  pool. If a `0_` s-exon has a protein sequence, Iduna includes and annotates it.
- A `0_` s-exon with no protein sequence adds no MSA columns, so it does not
  appear as a block in `*_s_exon_blocks.tsv`.

Each expansion PID directory also has `stage_state.json`, which records the
input identity used for safe cache reuse. If the seed MSA, expansion parameters,
or centroid request change, Iduna treats the cached expansion as stale, warns,
and rebuilds that PID directory. Validation directories use the same manifest
format and are reused when the seed, expansion, target, and UniProt comparison
inputs are unchanged.
