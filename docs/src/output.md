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
      centroid_msa/  # optional, written only with centroids=true / --centroids
  validation/
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

The final MSA paths are also available from the returned [`IdunaResult`](@ref):

```julia
result.thoraxe_msa.best_seed.stockholm_path
result.expansion.match_stockholm
result.expansion.a3m_path
```

`expanded_msa/` is the full MMseqs2 cluster-expanded MSA and remains the main
Iduna result. When `centroids=true` or `--centroids` is used, Iduna also writes
`centroid_msa/` with the MSA built from centroid/consensus hits before cluster
expansion. This side output reuses the same MMseqs2 search result, does not run
a second search, and does not change validation or `result.json`. Its files are
`<transcript>_centroids_full.sto`, `<transcript>_centroids_matchonly.sto`,
`<transcript>_centroids.a3m`, and `<transcript>_centroid_hits_raw.fasta`.
