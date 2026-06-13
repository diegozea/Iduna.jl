# Iduna 🍎

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/Iduna.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://diegozea.github.io/Iduna.jl/dev/)
[![Build Status](https://github.com/diegozea/Iduna.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/diegozea/Iduna.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/diegozea/Iduna.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/diegozea/Iduna.jl)
[![Coverage](https://coveralls.io/repos/github/diegozea/Iduna.jl/badge.svg?branch=main)](https://coveralls.io/github/diegozea/Iduna.jl?branch=main)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## What Iduna Does

Iduna _(Intrinsically Disordered Unit Aligner)_ builds a
[ThorAxe](https://github.com/diegozea/ThorAxe.jl) multiple sequence alignment
(MSA) for one UniProt accession or Ensembl transcript ID.

An MSA lines up related protein sequences so they can be compared residue by
residue. By default, Iduna expands the ThorAxe seed alignment with MMseqs2 and
HMMER, then writes validation summaries for the final result.

For transcript inputs, Iduna resolves the parent Ensembl gene and species before
calling ThorAxe.

## Quick Start

```julia
using Iduna

result = iduna("P20963"; mmseqs_db="/path/to/mmseqs/db")
expanded = load_expanded_msa(result)
```

`mmseqs_db` is required for the default full run. Use `no_expansion=true` when
you only need the ThorAxe seed MSA.

## Command Line Use

```bash
julia --project=. -m Iduna P20963 --mmseqs-db /path/to/mmseqs/db --workdir P20963
```

The installed `iduna` app accepts the same Iduna options. Julia thread flags
must be passed to Julia before the Iduna arguments; `--mmseqs-threads` controls
only the MMseqs2 expansion step.

## Main Outputs

Iduna writes a stable work directory. The main files and folders are:

- `result.json`: a summary that can be loaded later with `load_result`.
- `thoraxe_msa/`: ThorAxe seed alignments and seed-selection summaries.
- `expansion/`: expanded MSAs from MMseqs2/HMMER, unless expansion is skipped.
- `validation/`: size, diversity, and query-sequence checks.
- `logs/`: logs from the external tools.

Rerunning with `overwrite=false` reuses matching completed stages when their
inputs and required outputs still match.

## Common Options

### Stop After ThorAxe

Use `no_expansion=true` in Julia or `--no-expansion` in the app. In this mode,
`mmseqs_db` is not required.

```julia
thoraxe_only = iduna("ENST00000362089.10"; no_expansion=true,
    workdir="ENST00000362089_thoraxe")
seed = load_seed_msa(thoraxe_only)
```

### Choose Species

Iduna starts from `specieslist="ases"` by default, the curated species set used
by the Ases webserver. You can also pass `specieslist="all"`, a file path, a
comma-separated list, or one species name.

Smaller curated species lists can make the ThorAxe input step faster.

### Select Seeds

Iduna tests several ThorAxe percent identity (PID) thresholds and selects a seed
alignment. Set `pid_sample_count=0` to keep every eligible PID candidate instead
of selecting only one.

### Save Extra Files

Set `centroids=true` in Julia or `--centroids` in the app to also save a
centroid-level MSA. This is an extra output; the regular expanded MSA remains
the main result.

If a complete ThorAxe `transcript_query` bundle is already available, pass it as
`thoraxe_input_dir`; Iduna copies it into the work directory and still runs the
ThorAxe MSA stage.

## Documentation

Read the [stable documentation](https://diegozea.github.io/Iduna.jl/stable/) or
the [development documentation](https://diegozea.github.io/Iduna.jl/dev/) for
the full API, output layout, species filtering, seed selection, and reuse
behavior.
