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
           Base.find_package("ClustalO_jll") !== nothing
end
