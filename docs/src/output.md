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
  thoraxe/
    path_table.csv
    s_exon_table.csv
    msa/
  thoraxe_msa/
    seeds/
    best_seed.csv
  expansion/
    <gene>/<transcript>/
      seeds/
      dbs/
      hmm/
      logs/
      expanded_msa/
  validation/
    stats.csv
    query_vs_uniprot_alignment.txt
```

The output layout is part of the API. `overwrite=true` rebuilds package-owned
subdirectories, but Iduna does not delete files outside the selected work
directory.

When `thoraxe_input_dir` is supplied, that source bundle is treated as read-only
and copied into `thoraxe_input/`; the copied layout is then preserved like a
fresh `transcript_query` result.

The final MSA paths are also available from the returned [`IdunaResult`](@ref):

```julia
result.thoraxe_msa.best_seed.stockholm_path
result.expansion.match_stockholm
result.expansion.a3m_path
```
