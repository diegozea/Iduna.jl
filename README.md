# Iduna 🍎

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/Iduna.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://diegozea.github.io/Iduna.jl/dev/)
[![Build Status](https://github.com/diegozea/Iduna.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/diegozea/Iduna.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/diegozea/Iduna.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/diegozea/Iduna.jl)
[![Coverage](https://coveralls.io/repos/github/diegozea/Iduna.jl/badge.svg?branch=main)](https://coveralls.io/github/diegozea/Iduna.jl?branch=main)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Iduna _(Intrinsically Disordered Unit Aligner)_ builds a ThorAxe MSA for one UniProt 
accession or Ensembl transcript ID and expands the selected seed MSA with MMseqs2/HMMER. 
For transcript inputs, Iduna resolves the parent Ensembl gene and species before 
calling [ThorAxe](https://github.com/diegozea/ThorAxe.jl).

```julia
using Iduna

result = iduna("P20963"; mmseqs_db="/path/to/mmseqs/db", workdir="P20963")
expanded = load_expanded_msa(result)
```

The package writes a stable work directory containing ThorAxe outputs, PID seed
MSAs, expansion outputs, logs, and validation stats.

Pass `centroids=true` in Julia, or `--centroids` in the app, to also save a
centroid-level MSA before MMseqs2 expands centroid hits to cluster members. This
is a side output; the regular `expanded_msa/` files remain the main result.

By default, Iduna filters the ThorAxe species list with Ensembl homology using
`orthology="1:1"`, then checks BioMart Ensembl Gene dataset availability with
`biomart_datasets_filter=true`. BioMart dataset names are only used internally;
ThorAxe still receives species names. Use `orthology="1:n"` or `"m:n"` to keep
broader ortholog relationships, set `specieslist_filter=false` to skip the
Ensembl step, or set `biomart_datasets_filter=false` to skip the BioMart
dataset preflight.
`transcript_query_timeout_seconds` bounds each Ensembl download attempt and
Iduna can retry without the species list after a timeout. Set
`thoraxe_timeout_seconds` when individual ThorAxe runs should also have a
wall-clock limit.

If a complete ThorAxe `transcript_query` bundle is already available, pass it as
`thoraxe_input_dir`; Iduna copies it into the work directory and still runs the
ThorAxe MSA and MMseqs/HMMER expansion stages normally.
