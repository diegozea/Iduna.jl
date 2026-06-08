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
Iduna
iduna
load_result
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

## Internal Utility API

This section is intended for Iduna developers. These helpers support the public
pipeline API, tests, and maintenance tasks; they are not the main user-facing
API and may change as the internals evolve.

```@docs
Utils
Utils.DEFAULT_PID_THRESHOLDS
Utils.S_EXON_CODE_FEATURE
Utils.S_EXON_CODE_MAP_FEATURE
Utils.S_EXON_MISSING_CODE
Utils.strip_ensembl_version
Utils.sequence_name_variants
Utils.resolve_sequence_name
Utils.protein_alignment_stats
Utils.is_ensembl_transcript_id
Utils.is_uniprot_id
Utils.id_kind
Utils.format_pid
Utils.format_pid_dir
Utils.decode_body
Utils.prepare_output_dir
Utils.safe_rm
Utils.run_logged
Utils.fasta_sequence
Utils.format_fasta
Utils.write_fasta
Utils.write_text
Utils.write_json
Utils._pipeline_stage_dir
Utils._stage_state_path
Utils._write_stage_state
Utils._classify_stage_state
Utils.s_exon_blocks_path
Utils.s_exon_codes
Utils.has_s_exon_annotations
Utils.s_exon_code_map
Utils.set_s_exon_annotations!
Utils.write_s_exon_blocks_tsv
Utils.ensure_mmseqs_db
Utils.collect_stage_summaries
Utils.result_summary
```
