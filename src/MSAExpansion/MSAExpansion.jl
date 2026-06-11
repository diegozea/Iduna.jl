"""
    MSAExpansion

Expand a ThorAxe seed MSA with MMseqs2 and HMMER, then write expanded alignment
files for validation.
"""
module MSAExpansion

import HMMER_jll
import MMseqs2_jll

using Dates: UTC, now
using MIToS.MSA: A3M, AbstractMultipleSequenceAlignment, FASTASequences, Stockholm,
                 getannotcolumn, ncolumns, nsequences, read_file, sequence_id,
                 sequencenames, write_file

using ..Utils: ExpansionResult, ResolvedTarget, SeedSelection, _resolve_artifact_path,
               ensure_mmseqs_db, format_pid, format_pid_dir, run_logged, safe_rm,
               s_exon_code_map, s_exon_codes, set_s_exon_annotations!, write_fasta,
               write_s_exon_blocks_tsv, _classify_stage_state, _file_sha256,
               _read_stage_state, _stage_state_path, _write_stage_state

export expand_msa,
       normalize_stockholm_annotations!,
       prepare_stockholm_for_mmseqs

include("OutputPathsIdentity.jl")
include("CacheState.jl")
include("StockholmPreparation.jl")
include("CommandExecution.jl")
include("AlignmentProjection.jl")
include("HitCentroidOutputs.jl")
include("ExpansionWorkflow.jl")
include("PublicAPI.jl")

end
