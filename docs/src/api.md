# API

```@meta
CurrentModule = Iduna
```

## Main Entry Point

`iduna` runs the full ThorAxe plus MMseqs2/HMMER pipeline by default. Use
`no_expansion=true` to stop after the ThorAxe MSA stage; in that case
`isempty(result.expansions)`, `load_seed_msa(result)` remains available, and
`load_expanded_msa(result)` throws a guidance error.

```@docs
iduna
load_seed_msa
load_expanded_msa
```

## Result Types

```@docs
ResolvedTarget
SeedSelection
ThorAxeMSAResult
ExpansionResult
ValidationResult
IdunaResult
```

## Submodules

```@autodocs
Modules = [
    Iduna.IDMapping,
    Iduna.ThorAxeMSA,
    Iduna.MSAExpansion,
    Iduna.ResultsValidation,
]
```
