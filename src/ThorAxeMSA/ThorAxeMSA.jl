"""
    ThorAxeMSA

Build ThorAxe-based seed MSAs and choose the percent identity (PID) threshold
used by later Iduna stages.
"""
module ThorAxeMSA

import CSV
import Dates
import JSON
import ProgressMeter
import Scratch
import SHA
import ThorAxe

using Base.Threads
using DataFrames: DataFrame, nrow
using MIToS.MSA: AbstractMultipleSequenceAlignment, FASTA, PIRSequences, Stockholm,
                 getannotsequence, join_msas, nsequences, read_file, sequence_id,
                 sequencenames, setreference!, stringsequence, write_file
using Random: MersenneTwister
using StatsBase: sample

using ..EPLI
using ..Utils: DEFAULT_PID_THRESHOLDS, ResolvedTarget, SeedSelection, ThorAxeMSAResult,
               _http_get_request, _http_get_with_retries, _relative_artifact_path,
               _resolve_artifact_path, _sample_indices, _sample_rng,
               _terminal_progress_enabled,
               decode_body,
               fasta_sequence, format_pid, format_pid_dir, protein_alignment_stats,
               resolve_sequence_name, safe_rm,
               has_s_exon_annotations, s_exon_blocks_path, set_s_exon_annotations!,
               strip_ensembl_version, write_fasta, write_json, write_s_exon_blocks_tsv,
               write_text, _classify_stage_state, _pipeline_stage_dir,
               _read_stage_state, _write_stage_state

export assemble_transcript_msa,
       build_thoraxe_msa,
       compute_identity_against_reference,
       select_best_seed

include("ConstantsPaths.jl")
include("InputStageState.jl")
include("CommandExecution.jl")
include("SpeciesFiltering.jl")
include("TranscriptQuery.jl")
include("TranscriptMSAAssembly.jl")
include("TranscriptValidation.jl")
include("PIDCandidateArtifacts.jl")
include("CandidateScoring.jl")
include("CandidateSummaryCache.jl")
include("SeedSelection.jl")
include("ThorAxeMSAStage.jl")

end
