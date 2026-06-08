```@meta
CurrentModule = Iduna
```

# Iduna

## Overview

Iduna builds one ThorAxe-based multiple sequence alignment (MSA). By default, it
then expands that alignment with MMseqs2 and HMMER. The public API is centered
on [`iduna`](@ref), which accepts one UniProt accession or Ensembl transcript ID
and writes a reproducible work directory.

The package is file-first. External tools write logs and intermediate files
under the chosen `workdir`, and the returned result stores paths plus the
resolved identifiers, selected seed, expansion outputs, validation statistics,
warnings, and status. `result.workdir` is absolute in memory; saved metadata
uses relative paths for artifacts under `workdir`.

For an Ensembl transcript input, Iduna resolves the parent Ensembl gene ID and
species needed by ThorAxe. It does not require UniProt mapping on that path.

Rerunning with `overwrite=false` reuses matching completed stages. This includes
copied or previously generated `transcript_query` bundles. Iduna rebuilds a
stage when it is missing, stale, failed, or incomplete.

## Julia API

```julia
using Iduna

result = iduna(
    "P20963";
    mmseqs_db="/path/to/mmseqs/uniref_db",
    workdir="P20963",
    overwrite=false,
    centroids=false,
)

expanded = load_expanded_msa(result)
```

`mmseqs_db` is required for the default full run. Set `no_expansion=true` when
only the ThorAxe MSA stage is needed.

## Command Line

The same entry point is available as a Julia 1.12 app:

```bash
julia --project=. -m Iduna P20963 --mmseqs-db /path/to/mmseqs/uniref_db
```

The installed `iduna` app accepts the same Iduna options.

## ThorAxe-Only Runs

Use `no_expansion=true` in Julia, or `--no-expansion` in the app, to stop after
the ThorAxe MSA stage. In that mode `mmseqs_db` is not required,
`isempty(result.expansions)`, and `load_seed_msa(result)` loads the selected
ThorAxe percent identity (PID) seed.

```julia
thoraxe_only = iduna(
    "ENST00000362089.10";
    no_expansion=true,
    workdir="ENST00000362089_thoraxe",
)

seed = load_seed_msa(thoraxe_only)
thoraxe_only.thoraxe_msa.baseline_stockholms[1]
thoraxe_only.thoraxe_msa.seeds[1].stockholm_path
```

```bash
julia --project=. -m Iduna ENST00000362089.10 --no-expansion
```

## Centroid-Level MSA

Add `--centroids` (or `centroids=true` in Julia) to also save a
centroid-level MSA. This is a side output built from MMseqs2 centroid or
consensus hits before cluster expansion; the regular expanded MSA remains the
main result used by validation.

If an earlier cached expansion lacks the requested centroid files, Iduna warns
and rebuilds that PID expansion.

## Choosing Species

Iduna starts from `specieslist="ases"` by default, the curated set used on the
Ases webserver, filters it with Ensembl homology using `orthology="1:1"`, then
applies `biomart_datasets_filter=true` as a second preflight against the current
BioMart Ensembl Gene dataset list. BioMart dataset names are used only
internally; species names are still passed to ThorAxe. Use `specieslist="all"`
or `specieslist=""` for unrestricted ThorAxe species selection, pass a
comma-separated species list, file path, or single species name for an explicit
selection, use `orthology="1:n"` or `"m:n"` for broader ortholog relationships,
set `specieslist_filter=false` to skip the Ensembl homology step, or set
`biomart_datasets_filter=false` to skip the BioMart dataset preflight. The
lowercase strings `"ases"` and `"all"` are reserved presets; use `./ases` or
`./all` for files with those names.
The BioMart dataset list is cached in package scratch space and refreshed when
used on a later calendar date. Iduna also reports species recorded in
`transcript_query` BioMart failure outputs when a run completes with partial
BioMart failures.
`transcript_query` runtime depends heavily on the number of species it has to
download, so runs can be faster when you provide a small curated `specieslist`.
Use `thoraxe_input_dir` to reuse an existing complete `transcript_query` bundle
instead of fetching it again.

## Seed Selection

Seed selection is per percent identity (PID) threshold. Iduna runs
`transcript_query` once and builds one full-species candidate `msa_0` at each
PID threshold. Each candidate is validated: indels versus UniProt exclude that
PID from selection, while substitutions are reported as warnings. By default,
`sampling_strategy=:common`
draws one shared set of species samples from the species common to all eligible
`msa_0` candidates; each PID then runs ThorAxe with the same sampled species
lists and is scored by HHsuite against its own full `msa_0`. Use
`sampling_strategy=:independent` for the previous PID-local sampling behavior,
or `sampling_strategy=:input` to sample from the effective input species list.
Iduna chooses the highest median identity, highest mean identity, largest
candidate `msa_0` species count, and finally the first PID in `pid_thresholds`
order. By default Iduna uses 45 samples, retains 80% of non-reference species
per sample, and records a random `pid_sample_seed` unless one is supplied.
Set `pid_sample_count=0` to skip seed selection and carry every eligible PID
candidate forward. In that mode `result.thoraxe_msa.seeds`,
`result.expansions`, and `result.validations` can contain multiple entries.
`result.expansions` is indexed like `result.thoraxe_msa.seeds` and can contain
`missing` for an unavailable expansion result; use `pid=` or `index=` with
`load_seed_msa` and `load_expanded_msa` to select one.

## Threads

The repeated ThorAxe runs used for PID sample scoring can run at the same time.
For an installed `iduna` app, set Julia threads with an environment variable:

```bash
JULIA_NUM_THREADS=4 iduna P20963 --mmseqs-db /path/to/mmseqs/uniref_db --mmseqs-threads 8
```

or pass Julia flags before the app separator:

```bash
iduna --threads=4 -- P20963 --mmseqs-db /path/to/mmseqs/uniref_db --mmseqs-threads 8
```

From a source checkout, start Julia with more than one Julia thread:

```bash
julia --threads 4 --project=. -m Iduna P20963 --mmseqs-db /path/to/mmseqs/uniref_db --mmseqs-threads 8
```

This affects the ThorAxe sampling step. The app option `--mmseqs-threads` is
different: it controls MMseqs2 during the expansion step. Do not put
`--threads` after the Iduna arguments; it is a Julia runtime flag, not an Iduna
option.

## Reusing ThorAxe Input

If the ThorAxe `transcript_query` bundle has already been created, pass it with
`thoraxe_input_dir`. Iduna copies that bundle into `workdir/thoraxe_input` and
continues with the same ThorAxe MSA and percent identity (PID) seed stages. The
copied bundle is fingerprinted in a `thoraxe_input` stage manifest and reused on later
`overwrite=false` reruns when the target, species list, orthology, and filter
options match. Unless `no_expansion=true`, Iduna also runs expansion and
expanded-MSA validation.

## Output Details

The work directory contains ThorAxe outputs, PID seed MSAs, expansion outputs,
logs, validation statistics, and stage manifests. See [Output Layout](@ref) for
the full layout, path rules, and cache reuse details.

```@index
```
