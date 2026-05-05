```@meta
CurrentModule = Iduna
```

# Iduna

Iduna builds one ThorAxe-based multiple sequence alignment and expands it with
MMseqs2/HMMER. The public API is centered on [`iduna`](@ref), which accepts one
UniProt accession or Ensembl transcript ID and writes a reproducible work
directory.

The package is file-first. External tools write logs and intermediate files
under the chosen `workdir`, and the returned result stores paths plus the
resolved identifiers, selected seed, expansion outputs, validation statistics,
warnings, and status.

```julia
using Iduna

result = iduna(
    "P20963";
    mmseqs_db="/path/to/mmseqs/uniref_db",
    workdir="P20963",
    overwrite=false,
    transcript_query_timeout_seconds=180,
)

expanded = load_expanded_msa(result)
```

The same entry point is available as a Julia 1.12 app:

```bash
julia --project=. -m Iduna P20963 --mmseqs-db /path/to/mmseqs/uniref_db
```

For an Ensembl transcript input, Iduna resolves the parent Ensembl gene ID and
species needed by ThorAxe. It does not require UniProt mapping on that path.

`transcript_query_timeout_seconds` defaults to 180 seconds, with a bounded retry
that can drop the species list after a timeout. `thoraxe_timeout_seconds` is
unset by default because ThorAxe runtime depends on gene complexity and the
selected PID thresholds.

If the ThorAxe `transcript_query` bundle has already been created, pass it with
`thoraxe_input_dir`. Iduna copies that bundle into `workdir/thoraxe_input` and
continues with the same ThorAxe MSA, PID seed, expansion, and validation stages.

```@index
```
