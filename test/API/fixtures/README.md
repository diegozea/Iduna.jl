# Legacy custom-PID archive

`pre59-custom-pid/` is an unchanged output archive written by the predecessor of
the precise PID naming change, commit `b067e488796644b4206bc93551f7f6a4ba6aa680`.
It used the ordinary result writer and validation stage with deterministic target
and ThorAxe callbacks. The selected PID is `30.001`, while its saved artifacts use
the historical `pid_30.00` directory. It intentionally has no candidate-summary CSV.

The `result.json` SHA-256 is
`8b83d4409ee23b594794f30d08ad67c302dc061b836afe96022ebbba2b77e350`.
Keep all 11 archive files unchanged so tests remain independent of the current
writer and path formatter. Loader tests also cover relocation and coexistence with
the newer precise directory spelling.
