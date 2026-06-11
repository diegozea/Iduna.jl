"""
    EPLI

Estimate profile-level identity from repeated MSA subsamples.
"""
module EPLI

import CSV
import HHsuite_jll
import JSON
import ProgressMeter
import ThorAxe

using Base.Threads
using DataFrames: DataFrame
using MIToS.MSA: FASTA, FASTASequences, nsequences, read_file, sequence_id,
                 stringsequence
using Statistics: median

using ..Utils: _file_sha256, _sample_indices, _sample_rng, _terminal_progress_enabled,
               write_fasta, write_json

export comparable_positions_normalization,
       epli_score,
       hhsuite_identity_score,
       hhsuite_profile_score,
       no_normalization,
       self_reference_normalization

include("HHsuite.jl")
include("Normalization.jl")
include("AlignmentCache.jl")
include("Sampling.jl")
include("Scoring.jl")
include("JLLAligners.jl")
include("ProGraphMSA.jl")

end
