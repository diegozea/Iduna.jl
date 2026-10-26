# EPLI Score

EPLI means Estimated Profile-Level Identities. It measures how stable an aligner
is when the input sequence set is subsampled.

## Quick Example

This example creates four fake protein sequences with MIToS, aligns them with
ProGraphMSA, and computes EPLI with the default sampling settings.

```jldoctest; output = false
using Iduna
using MIToS.MSA

sequences = [
    AnnotatedSequence("ref", "MKTAYIAKQRQISFVKSHFSRQDILD"),
    AnnotatedSequence("seq2", "MKTAYIAKQRQISFVKSHFSRQNILD"),
    AnnotatedSequence("seq3", "MKTAYIAKQKQISFVKSHFSRQDILD"),
    AnnotatedSequence("seq4", "MKTAYIAKQRQISFVKSHFTRQDILD"),
]

workdir = mktempdir()

result = Iduna.EPLI.epli_score(
    sequences,
    workdir,
    Iduna.EPLI.prographmsa_aligner,
)

(result.n_samples, result.input_source, isfile(result.summary_path))

# output

[ Info: EPLI alignment reuse disabled: no aligner_cache_id was supplied.
(45, "mitos_sequences", true)
```

The call writes a full input FASTA under `workdir`, runs ProGraphMSA for the full
sequence set and for each sampled subset, compares each sampled MSA to the full
MSA with HHsuite, and stores the median normalized score in `result.epli`.
This temporary-directory example does not enable persistent alignment reuse, so
Iduna reports an informational cache notice.

The main output files are:

- `result.summary_path`: one-row CSV with the run settings and final EPLI.
- `result.scores_path`: per-sample scores used to compute the median.
- `result.reference_msa_fasta`: the full ProGraphMSA alignment.

## Input Types

`epli_score` accepts these sequence sources:

- An unaligned FASTA path.
- A `Vector` of MIToS `AbstractSequence` objects, such as `AnnotatedSequence`.
- A single MIToS `AbstractSequence`.
- A `MIToS.MSA.AbstractMultipleSequenceAlignment`.

MIToS sequence and MSA inputs are converted to unaligned FASTA before the aligner
runs. Gap characters `-` and `.` are removed from each sequence.

For aligned MSA input:

```julia
using Iduna
using MIToS.MSA

msa = read_file("aligned_sequences.fasta", FASTA)

result = Iduna.EPLI.epli_score(
    msa,
    "msa_epli",
    Iduna.EPLI.prographmsa_aligner;
    reference_sequence = "seq1",
)
```

The returned result and `summary.csv` use `input_fasta = missing` for MIToS
object inputs.

## Wrap An Aligner

An aligner wrapper is a Julia function with this required contract:

```julia
function my_aligner(input_fasta, output_fasta; logs_dir = nothing,
        run_label = "run", aligner_args = Cmd(String[]))
    run(`my-aligner --input $input_fasta --output $output_fasta $aligner_args`)
    return output_fasta
end
```

The wrapper receives an unaligned FASTA and must write an aligned FASTA. It can
return the output path, return `nothing`, or return a named tuple with a
`fasta_path` field.

EPLI sets these arguments when it calls the wrapper:

- `input_fasta`: unaligned FASTA for this run. It can contain the full input set
  or one sampled subset.
- `output_fasta`: path where EPLI expects the aligned FASTA for this run.
- `logs_dir`: optional directory for logs.
- `run_label`: stable label for the current run, such as `full` or
  `sequence_subset_001`.
- `aligner_args`: extra command arguments supplied by the `epli_score` caller.

`aligner_args` is intentionally a `Cmd`, not open-ended Julia keywords. This
keeps EPLI options separate from command-line options for the external aligner.
Wrappers must accept `run_label` and `aligner_args`, even when callers use the
default empty command.

## Reusing Alignments

By default, `epli_score` runs the aligner again. To enable persistent alignment
reuse, supply a nonblank `aligner_cache_id` that identifies the callback code,
external tool build, optional dependencies and any other state that affects its
output. Change this ID when any of those inputs changes. Iduna trusts this caller
assertion; a label or function name alone does not establish provenance. The same
rule applies to the bundled aligner wrappers.

Reuse also requires known Iduna source and Julia dependency identities, matching
input sequences, arguments and run labels, and verified output contents. The
implementation record is saved in each alignment's `.iduna_cache.json` sidecar.
An unavailable identity disables reuse; computation remains available. Old
function-hash caches rebuild once a trustworthy identity is supplied. Scores are
always recomputed, and sampling and the EPLI formula are unchanged.

Source files and loaded dependencies must stay unchanged for the Julia process;
restart Julia after editing them. A fixed `sample_seed` reproduces sequence-row
subsets, but does not make a stochastic external aligner deterministic. Include
relevant tool seeds and configuration in the callback identity or its arguments.

## Optional JLL Aligners

Iduna also provides unexported wrappers for FAMSA, Kalign, MAFFT, Clustal
Omega, and MUSCLE through Julia package extensions. These wrappers are
available only after loading the matching JLL package.

All bundled aligner wrappers return aligned FASTA records in the same order as
the input FASTA. This keeps EPLI's reference row stable across aligners.
Wrappers that reorder aligner output use MIToS FASTA parsing and writing, so
input and aligned output sequences should use the MIToS residue alphabet.

For MAFFT:

```julia
using Iduna
using MAFFT_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "mafft_epli",
    Iduna.EPLI.mafft_aligner,
    aligner_args = `--thread 4 --localpair --maxiterate 1000`,
)
```

For FAMSA:

```julia
using Iduna
using FAMSA_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "famsa_epli",
    Iduna.EPLI.famsa_aligner,
    aligner_args = `-t 4`,
)
```

For Kalign:

```julia
using Iduna
using kalign_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "kalign_epli",
    Iduna.EPLI.kalign_aligner,
    aligner_args = `-n 4`,
)
```

The Kalign wrapper fixes FASTA output with `--format fasta`. It does not force
`--type protein`, because Kalign rejects DNA-like protein alphabets such as
`AAAA` before alignment when that flag is set. Pass `--type` in `aligner_args`
only when you want Kalign to validate a specific sequence type.

For Clustal Omega:

```julia
using Iduna
using ClustalO_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "clustalo_epli",
    Iduna.EPLI.clustalo_aligner,
    aligner_args = `--threads=4 --auto`,
)
```

For MUSCLE:

```julia
using Iduna
using MUSCLE_jll

result = Iduna.EPLI.epli_score(
    "sequences.fasta",
    "muscle_epli",
    Iduna.EPLI.muscle_aligner,
)
```

Use the qualified names `Iduna.EPLI.famsa_aligner`,
`Iduna.EPLI.kalign_aligner`, `Iduna.EPLI.mafft_aligner`,
`Iduna.EPLI.clustalo_aligner`, and `Iduna.EPLI.muscle_aligner`; they are
intentionally not exported.

## Alternative Scores

By default, EPLI uses HHsuite identity counts and comparable-position
normalization:

```text
100 * matched_positions / comparable_positions
```

For profile-score experiments, pass a different score and normalization function:

```julia
using Iduna
using Iduna.EPLI

result = epli_score(
    "sequences.fasta",
    "epli_work",
    Iduna.EPLI.prographmsa_aligner;
    score_fn = hhsuite_profile_score,
    normalization_fn = self_reference_normalization,
)
```

This computes each sample's raw HHalign profile score and normalizes it by the
raw score from aligning the full reference MSA against itself.

Custom score functions must return a finite `raw_score`. Normalization functions
must return a finite score, or a named tuple with a finite `normalized_score` and
an optional finite or `missing` `normalization_score`. These values must remain
finite when converted to `Float64`. Invalid scores raise `ArgumentError`
(wrapped in a task failure during parallel scoring) rather than removing samples
from the median. Finite negative scores and ratios above 100 are allowed; self
normalization still rejects a zero reference score.
Seed selection also rejects non-finite scores in eligible or cached selected
candidates. Missing scores remain valid when PID sampling is explicitly disabled.

Iduna uses EPLI internally to choose the best seed generated by ThorAxe. That
adapter uses the same scoring layer, but its sampling step is different: ThorAxe
does not normally accept arbitrary sampled sequence rows as input, so Iduna
samples species lists instead of sequences for ThorAxe PID selection.

```@docs
Iduna.EPLI.epli_score
Iduna.EPLI.hhsuite_identity_score
Iduna.EPLI.hhsuite_profile_score
Iduna.EPLI.comparable_positions_normalization
Iduna.EPLI.self_reference_normalization
Iduna.EPLI.no_normalization
Iduna.EPLI.prographmsa_aligner
```
