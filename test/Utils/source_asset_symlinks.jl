using Iduna, Test

@testset "Symlinked scientific assets" begin
    mktempdir() do root
        fingerprint = Iduna.Utils._source_fingerprint
        write(joinpath(root, "Project.toml"), "name = \"LinkedAssets\"\n")
        data = mkpath(joinpath(root, "scientific_data"))
        matrix = joinpath(data, "matrix")
        write(matrix, "A 4\n")
        for directory in ("src", "ext")
            tree = mkpath(joinpath(root, directory))
            for name in ("data", "Tiny.jl.cov", "Tiny.jl.123.cov")
                linked_directory = joinpath(tree, name)
                symlink(data, linked_directory)
                @test fingerprint(root) === nothing
                rm(linked_directory)
            end
            for name in ("matrix", "matrix.cov")
                linked_file = joinpath(tree, name)
                symlink(matrix, linked_file)
                before = fingerprint(root)
                @test before !== nothing
                write(matrix, "A -99\n")
                @test fingerprint(root) != before
                write(matrix, "A 4\n")
                @test fingerprint(root) == before
                rm(linked_file)
            end
        end
    end
end
