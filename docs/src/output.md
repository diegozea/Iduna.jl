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
directory. Cleanup checks both the path and its resolved parent, and refuses to
remove the work directory itself. It can remove an in-directory symlink without
deleting its target, but rejects cleanup through a parent symlink that leaves
the work directory. A symlinked work directory uses its resolved target as the
boundary. These checks assume directories and symlinks are not changed concurrently.

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

Completed manifests also record SHA-256 digests for every required scientific file.
Reuse requires the current schema, matching inputs and output paths, and unchanged
contents of all those files. Empty raw-hit files are valid; final alignments must
parse and contain sequences. Logs and reproducible s-exon block tables are excluded
from the required digest set; missing block tables are rebuilt from verified alignments.
A result directory can be moved without invalidating its internal file paths.

Older manifests, metadata-only caches and partial outputs are not adopted. The next
requested run rebuilds the affected package-owned stage and may replace its previous
outputs. There is no automatic archive migration or global cache deletion. Untouched
archived copies remain available through `load_result`, which reads historical files
without certifying their provenance. Use one writer per work directory; these checks
do not synchronize simultaneous runs or other programs modifying the same files.

Custom EPLI aligners must return a nonempty aligned FASTA with unique full record IDs
matching the input, equal row lengths, and unchanged ungapped residues. Case and gap
characters (`-` or `.`) do not affect the residue comparison; row reordering, wrapping,
and alternate returned file paths are supported. A failed callback leaves no reusable
completion marker for its partially written output.

`thoraxe_msa/runs/` keeps the raw ThorAxe output for each percent identity (PID)
candidate and its species samples. `thoraxe_msa/candidates/` stores each PID's
full-species candidate `msa_0`, sampled species list paths, gap-free sequence
files, and per-PID score CSV. With `sampling_strategy=:common` or `:input`,
shared sampled species lists live in `thoraxe_msa/samples/species/` and each
PID-local sampled species path is a symlink to the matching shared file.
`thoraxe_msa/candidate_summary.csv` records validation status, eligibility,
EPLI scores, candidate size, and the selected seed rows.

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

By default, Iduna selects one percent identity (PID) seed. Set
`pid_sample_count=0` to skip seed selection and carry every eligible PID
candidate forward, producing one
validation directory and, unless `no_expansion=true`, one expansion directory
per selected PID. `result.expansions` is seed-indexed and may contain `missing`
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
| `pid` | the ThorAxe percent identity (PID) threshold used for that seed |
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

Each expansion percent identity (PID) directory also has `stage_state.json`,
which records the input identity used for safe cache reuse. If the seed MSA,
expansion parameters, or centroid request change, Iduna treats the cached
expansion as stale, warns, and rebuilds that PID directory. Validation
directories use the same manifest format and are reused when the seed,
expansion, target, and UniProt comparison inputs are unchanged.


## Implementation and tool provenance

Reuse records the Iduna source fingerprint, Julia version, selected and loaded Julia
dependency versions and source trees, and relevant external tool identities. Scientific
stages include the hard dependency closure, rather than only direct package versions.
Path-developed packages hash their `Project.toml` and all files under `src` and `ext`,
including scientific data such as substitution matrices. Generated Julia coverage
files (`*.jl.cov` and `*.jl.<pid>.cov`) are excluded.
File symlinks hash their target contents. Untraversed directory symlinks make the
source identity unknown and disable affected reuse.
Registered package and JLL artifact trees are assumed immutable. Unknown package,
source or executable identity disables the affected reuse, while computation remains
available. Restart Julia after source or dependency changes; hot reload and concurrent
source editing during a run are outside this contract.

ThorAxe provenance includes its actual Python interpreter, installed ThorAxe source
and distribution dependencies, CLI scripts, ProGraphMSA, and HHsuite where used.
The Python probe retains the CLI's interpreter invocation path and records its
environment prefixes, including a virtual environment reached through a symlink.
Different interpreter paths in the two CLI shebangs disable affected reuse.
Python dependency versions and locations identify installed trees; their contents
and the Python installation must remain immutable for those recorded versions.
In-place edits to dependency code are not detected by version metadata.
A copied input bundle is identified by its contents and does not require local Python
provenance merely to reuse that bundle. MSA selection also binds reference sequence
contents and species-list file contents. Generated sampling seeds retain their recorded
value; other missing metadata is not treated as a match. CPU, Julia thread count,
relevant thread environment variables, and explicit MMseqs thread settings participate
in reuse. An unspecified MMseqs thread count cannot authorize expansion reuse.

Generic EPLI callbacks require an explicit `aligner_cache_id` for persistent alignment
reuse, including built-in callbacks. The caller must identify the callback implementation,
external tool build and relevant external state. Without an ID each call realigns.
See [EPLI](epli.md) for the callback contract; the EPLI formula is unchanged.

Legacy completion records and candidate caches require recomputation, including those
whose producer fingerprints predate these interpreter and source-data checks. Iduna does not
rename archives or delete all caches as a migration, but rerunning a stale stage replaces
its package-owned outputs. Copy results you need to retain before rerunning. Loading
recorded paths from untouched archives does not certify their historical scores.

## Immutable MMseqs database snapshots

Expansion reuse requires `mmseqs_db_manifest="/path/to/snapshot.json"` in Julia or
`--mmseqs-db-manifest /path/to/snapshot.json` in the app. The manifest is a JSON object
with `schema_version: 1`, `immutable: true`, a nonempty `snapshot_id`, and `components`
mapping every database component filename to its lowercase SHA-256 checksum.
Components include the representative, `_seq`, `_aln`, `_h` and `_seq_h` families,
their data, `.index` and `.dbtype` files, and any other files belonging to that prefix.
Split data starts at `.0`. The manifest itself is excluded from its inventory.

Prepare the manifest once, after database creation has finished, and keep the database
immutable for as long as its results may be reused. Iduna checks the declaration and
current file inventory and binds the complete manifest to the expansion identity. It
does **not** rehash a large database on every call. An unchanged manifest is not proof
that someone has left the underlying bytes unchanged: that is the operator's explicit
immutability obligation. After changing any component, prepare a new manifest before
reuse. With no valid manifest, expansion computes again instead of reusing saved results.

The inventory consists of files named exactly like the database prefix or beginning
with that basename followed by `.` or `_`. Keep unrelated logs elsewhere. If components
are symlinks, both links and targets must remain immutable. The following standalone
Python example streams each file once; change its three configuration values first.
It uses exclusive creation to avoid replacing an existing snapshot declaration.

```python
import hashlib
import json
from pathlib import Path

prefix = Path("/data/mmseqs/uniref_release/uniref").absolute()
manifest = Path("/data/manifests/uniref_release.json").absolute()
snapshot_id = "uniref-release-build-1"

components = {}
for path in sorted(prefix.parent.iterdir()):
    matches = (path.name == prefix.name
               or path.name.startswith(prefix.name + ".")
               or path.name.startswith(prefix.name + "_"))
    if not matches or not path.is_file() or path == manifest:
        continue
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    components[path.name] = digest.hexdigest()

assert components, "No database components found"
manifest.parent.mkdir(parents=True, exist_ok=True)
# Exclusive creation avoids replacing an existing snapshot declaration.
with manifest.open("x") as stream:
    json.dump({"schema_version": 1, "immutable": True,
               "snapshot_id": snapshot_id, "components": components},
              stream, indent=2, sort_keys=True)
    stream.write("\n")
```

Pass the manifest with its database prefix:

```julia
iduna("P28482";
    workdir = "/path/to/work",
    mmseqs_db = "/data/mmseqs/uniref_release/uniref",
    mmseqs_db_manifest = "/data/manifests/uniref_release.json",
    mmseqs_threads = 4)
```
