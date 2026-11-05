function write_epli_test_input(dir; fasta = ">ref\nACDE\n>other\nACDF\n",
        filename = "sequences.fasta")
    path = joinpath(dir, filename)
    write(path, fasta)
    return path
end

# Construct inside the callback so invalid values fail at the scoring boundary.
constant_score_callback(value) = (args...; kwargs...) -> Iduna.EPLI.RawScore(value)

# These test types deliberately accept unchecked data to exercise callback validation.
struct CustomRawScore{V} <: Iduna.EPLI.ScoreResult
    value::V
end
Iduna.EPLI.raw_score(score::CustomRawScore) = score.value

struct CustomCountScore{C, D} <: Iduna.EPLI.ScoreResult
    value::Float64
    evidence::C
    details::D
end
Iduna.EPLI.raw_score(score::CustomCountScore) = score.value
Iduna.EPLI.identity_counts(score::CustomCountScore) = score.evidence
Iduna.EPLI.score_diagnostics(score::CustomCountScore) = score.details

function sample_size_score(reference_msa, sample_msa; kwargs...)
    n_sample = Iduna.EPLI.nsequences(Iduna.EPLI.read_file(sample_msa, Iduna.EPLI.FASTA))
    return Iduna.EPLI.RawScore(n_sample)
end

function copy_aligner(input_fasta, output_fasta; logs_dir = nothing,
        run_label = "run", aligner_args = Cmd(String[]))
    mkpath(dirname(output_fasta))
    cp(input_fasta, output_fasta; force = true)
    return output_fasta
end

function padded_aligner(input_fasta, output_fasta; logs_dir = nothing,
        run_label = "run", aligner_args = Cmd(String[]))
    records = Tuple{String, String}[]
    current_name = nothing
    current_seq = String[]
    for line in eachline(input_fasta)
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, '>')
            if current_name !== nothing
                push!(records, (current_name, join(current_seq)))
            end
            current_name = stripped[2:end]
            empty!(current_seq)
        else
            push!(current_seq, stripped)
        end
    end
    if current_name !== nothing
        push!(records, (current_name, join(current_seq)))
    end
    width = maximum(length(seq) for (_name, seq) in records)
    mkpath(dirname(output_fasta))
    open(output_fasta, "w") do io
        for (name, seq) in records
            println(io, ">", name)
            println(io, rpad(seq, width, '-'))
        end
    end
    return output_fasta
end

function jll_extension_test_dependencies_available()
    return Base.find_package("MAFFT_jll") !== nothing &&
           Base.find_package("ClustalO_jll") !== nothing &&
           Base.find_package("FAMSA_jll") !== nothing &&
           Base.find_package("MUSCLE_jll") !== nothing &&
           Base.find_package("kalign_jll") !== nothing
end
