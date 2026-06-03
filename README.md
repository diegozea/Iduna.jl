# Iduna 🍎

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/Iduna.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://diegozea.github.io/Iduna.jl/dev/)
[![Build Status](https://github.com/diegozea/Iduna.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/diegozea/Iduna.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/diegozea/Iduna.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/diegozea/Iduna.jl)
[![Coverage](https://coveralls.io/repos/github/diegozea/Iduna.jl/badge.svg?branch=main)](https://coveralls.io/github/diegozea/Iduna.jl?branch=main)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Iduna _(Intrinsically Disordered Unit Aligner)_ builds a ThorAxe MSA for one UniProt
accession or Ensembl transcript ID and, by default, expands the selected seed MSA
with MMseqs2/HMMER.
For transcript inputs, Iduna resolves the parent Ensembl gene and species before
calling [ThorAxe](https://github.com/diegozea/ThorAxe.jl).

```julia
using Iduna

result = iduna("P20963"; mmseqs_db="/path/to/mmseqs/db", workdir="P20963")
expanded = load_expanded_msa(result)
```

The package writes a stable work directory containing ThorAxe outputs, PID seed
MSAs, expansion outputs, logs, and validation stats. `result.workdir` is
absolute, while artifact paths under it are reported relative to `workdir`.
Pipeline stages record `stage_state.json` manifests with input identities and
required outputs. Rerunning with `overwrite=false` reuses completed matching
stages, including `transcript_query` bundles, and rebuilds only missing, stale,
failed, or incomplete stages.

For seed selection, Iduna runs `transcript_query` once and builds one
full-species candidate `msa_0` at each PID threshold. Each candidate is validated:
indels versus UniProt exclude that PID from selection, while substitutions are
reported as warnings. By default, `sampling_strategy=:common` draws one shared
set of species samples from the species common to all eligible candidates, runs
each PID with those same sampled species lists, and scores by HHsuite against
that PID's own full `msa_0`. Use `sampling_strategy=:independent` for PID-local
sampling or `sampling_strategy=:input` to sample from the effective input species
list. The selected seed is the highest median identity, then highest mean
identity, then the largest candidate `msa_0` species count, then the first PID
in `pid_thresholds` order. The default is 45 samples, 80% of non-reference
species per sample, and a random `pid_sample_seed` recorded in `result.json`.

The repeated ThorAxe runs used for PID sample scoring can run at the same time.
Start Julia with more than one Julia thread before running Iduna:

```bash
julia --threads 4 --project=. -m Iduna P20963 --mmseqs-db /path/to/mmseqs/db --mmseqs-threads 8
```

You can also set `JULIA_NUM_THREADS=4` before starting Julia. This affects the
ThorAxe sampling step. The app option `--mmseqs-threads` is different: it
controls MMseqs2 during the expansion step.

To stop after the ThorAxe MSA stage, use `no_expansion=true` in Julia or
`--no-expansion` in the app. In that mode `mmseqs_db` is not required,
`isempty(result.expansions)`, and the ThorAxe MSA paths are available from
`result.thoraxe_msa.baseline_stockholms`,
`result.thoraxe_msa.baseline_fastas`,
`result.thoraxe_msa.seeds[1].stockholm_path`, and
`result.thoraxe_msa.seeds[1].fasta_path`.

```julia
thoraxe_only = iduna("ENST00000362089.10"; no_expansion=true,
    workdir="ENST00000362089_thoraxe")
seed = load_seed_msa(thoraxe_only)
```

Set `pid_sample_count=0` to skip PID seed selection and carry every eligible
PID candidate forward as a seed. In that mode `result.thoraxe_msa.seeds`,
`result.expansions`, and `result.validations` may contain more than one entry;
pass `pid=` or `index=` to `load_seed_msa` or `load_expanded_msa` when needed.

Pass `centroids=true` in Julia, or `--centroids` in the app, to also save a
centroid-level MSA before MMseqs2 expands centroid hits to cluster members. This
is a side output; the regular `expanded_msa/` files remain the main result.
If an earlier cached expansion lacks the requested centroid files, Iduna warns
and rebuilds that PID expansion.

By default, Iduna filters the ThorAxe species list with Ensembl homology using
`orthology="1:1"`, then checks BioMart Ensembl Gene dataset availability with
`biomart_datasets_filter=true`. BioMart dataset names are only used internally;
ThorAxe still receives species names. Use `orthology="1:n"` or `"m:n"` to keep
broader ortholog relationships, set `specieslist_filter=false` to skip the
Ensembl step, or set `biomart_datasets_filter=false` to skip the BioMart
dataset preflight.
`transcript_query` runtime depends heavily on the number of species it has to
download, so runs can be faster when you provide a small curated `specieslist`.

If a complete ThorAxe `transcript_query` bundle is already available, pass it as
`thoraxe_input_dir`; Iduna copies it into the work directory and still runs the
ThorAxe MSA stage. Unless `no_expansion=true`, it then continues with the
MMseqs/HMMER expansion stages normally. Copied or previously generated
`thoraxe_input/Ensembl` bundles are fingerprinted and reused on later reruns
when the target, species list, orthology, and filter options match.
