"""
    Utils

Shared helpers for IDs, files, alignment checks, stage summaries, and s-exon
annotations used by Iduna.
"""
module Utils

using BioAlignments: AffineGapScoreModel, BLOSUM62, GlobalAlignment, alignment,
                     count_deletions, count_insertions, count_mismatches, pairalign
using BioSequences: LongAA
import CodecZlib
using Dates: UTC, now
import HTTP
import JSON
using Random: MersenneTwister
import SHA
using MIToS.MSA: AbstractMultipleSequenceAlignment, getannotcolumn, getannotfile,
                 ncolumns, sequencenames, setannotcolumn!, setannotfile!
using Printf: @sprintf
using StatsBase: sample

export DEFAULT_PID_THRESHOLDS,
       ExpansionResult,
       IdunaResult,
       ResolvedTarget,
       SeedSelection,
       ThorAxeMSAResult,
       ValidationResult,
       decode_body,
       ensure_mmseqs_db,
       fasta_sequence,
       format_fasta,
       id_kind,
       is_ensembl_transcript_id,
       is_uniprot_id,
       prepare_output_dir,
       protein_alignment_stats,
       resolve_sequence_name,
       result_summary,
       run_logged,
       safe_rm,
       sequence_name_variants,
       strip_ensembl_version,
       write_fasta,
       write_json,
       write_text,
       format_pid,
       format_pid_dir,
       S_EXON_CODE_FEATURE,
       S_EXON_CODE_MAP_FEATURE,
       S_EXON_MISSING_CODE,
       s_exon_blocks_path,
       s_exon_code_map,
       s_exon_codes,
       has_s_exon_annotations,
       set_s_exon_annotations!,
       write_s_exon_blocks_tsv

include("Types.jl")
include("Constants.jl")
include("TerminalSampling.jl")
include("HTTP.jl")
include("Identifiers.jl")
include("Files.jl")
include("ArtifactPaths.jl")
include("StageState.jl")
include("SExonAnnotations.jl")
include("ExternalDatabases.jl")
include("ResultSummaries.jl")

end
