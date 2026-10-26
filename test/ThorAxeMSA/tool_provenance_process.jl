using Iduna, Test

# A separate process keeps wrapper references, including initially unassigned refs,
# isolated from the rest of the suite.
@testset "Installed ThorAxe producer identity" begin
    msa = Iduna.ThorAxeMSA
    wrapper = msa.ThorAxe
    # Use a fresh directory so cached tools cannot hide missing setup.
    wrapper.ALIGNER_PATH[] = mktempdir()
    wrapper.aligner_executable()
    installed_aligner_path = wrapper.ALIGNER_PATH[]
    tools = msa._thoraxe_tool_identity(; aligner = false)
    @test tools !== nothing
    if tools !== nothing
        metadata = tools.python.metadata
        @test occursin(r"^[0-9a-f]{64}$", metadata["source_sha256"])
        @test !isempty(metadata["import_paths"])
        @test haskey(metadata["dependency_versions"], "thoraxe")
        @test keys(metadata["dependency_versions"]) ==
              keys(metadata["dependency_locations"])
    end
    for pythonpath in (".", "", "/tmp:")
        withenv("PYTHONPATH" => pythonpath) do
            @test_logs (:warn, r"ThorAxe tool identity unavailable") begin
                @test msa._thoraxe_tool_identity(; aligner = false) === nothing
            end
        end
    end

    sampled = msa._thoraxe_stage_provenance(1)
    @test msa._known_thoraxe_provenance(sampled, 1)
    @test sampled.python == tools.python
    @test sampled.aligner == msa._thoraxe_aligner_identity()
    @test Set(keys(sampled.hhsuite.executables)) == Set(["hhmake", "hhalign"])
    unsampled = msa._thoraxe_stage_provenance(0)
    @test msa._known_thoraxe_provenance(unsampled, 0)
    @test unsampled.python == sampled.python
    @test unsampled.hhsuite === nothing

    mktempdir() do root
        original = joinpath(wrapper.ALIGNER_PATH[], wrapper.ALIGNER_NAME)
        @test isfile(original)
        copied = joinpath(mkpath(joinpath(root, "with spaces")), wrapper.ALIGNER_NAME)
        cp(original, copied)
        wrapper.ALIGNER_PATH[] = dirname(copied)
        wrapper.ALIGNER_SHIM_PATH[] = joinpath(root, "missing-launcher")
        @test_throws ErrorException msa._thoraxe_aligner_identity()

        launcher = joinpath(root, "launcher")
        symlink(copied, launcher)
        wrapper.ALIGNER_SHIM_PATH[] = launcher
        identity = msa._thoraxe_aligner_identity()
        @test identity.path == realpath(copied)
        @test identity.sha256 == bytes2hex(Iduna.Utils.SHA.sha256(read(original)))
        open(copied, "a") do io
            write(io, "changed")
        end
        @test_throws ErrorException msa._thoraxe_aligner_identity()
    end
    wrapper.ALIGNER_PATH[] = installed_aligner_path
end

include("venv_provenance.jl")
include("relative_path_provenance.jl")
