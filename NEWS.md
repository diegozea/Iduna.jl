## Iduna.jl Release Notes

### Changes from v0.1.0 to v0.2.0

- Added `no_expansion=true` and `--no-expansion` to stop after the ThorAxe MSA
  stage without requiring an MMseqs2 database; in that mode,
  `IdunaResult.expansion` is `nothing` and the JSON summary stores
  `"expansion": null`.
- Added `centroids=true` and `--centroids` to write the centroid-level MSA side
  output before MMseqs2 expands hits to all cluster members.
- Added one-field-per-line pretty printing for `IdunaResult` and the nested
  result objects returned by the pipeline.
