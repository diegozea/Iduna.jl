# Actual tool provenance
# ----------------------

const _THORAXE_IDENTITY_PROBE = raw"""
import hashlib, json, os, sys
from pathlib import Path
python_environment = {name: os.environ.get(name) for name in ('PYTHONPATH', 'PYTHONHOME')}
pythonpath = python_environment['PYTHONPATH']
if pythonpath is not None and any(not Path(p).is_absolute() for p in pythonpath.split(os.pathsep)):
    raise RuntimeError('Relative PYTHONPATH cannot identify tools across working directories.')
sys.path[0] = str(Path(sys.argv[1]).resolve().parent)
import thoraxe, pkg_resources
from thoraxe.version import __version__
root = Path(thoraxe.__file__).resolve().parent
files = sorted(p for p in root.rglob('*') if p.is_file()
               and '__pycache__' not in p.parts and p.suffix not in ('.pyc', '.pyo'))
records = ''.join(p.relative_to(root).as_posix() + '\0'
                  + hashlib.sha256(p.read_bytes()).hexdigest() + '\n' for p in files)
dist = pkg_resources.get_distribution('thoraxe')
dependencies = pkg_resources.working_set.resolve([dist.as_requirement()])
print(json.dumps(dict(
    python_executable=sys.executable, python_version=sys.version,
    python_prefix=sys.prefix, python_base_prefix=sys.base_prefix,
    import_paths=sys.path, python_environment=python_environment,
    thoraxe_version=__version__, distribution_version=dist.version,
    source_sha256=hashlib.sha256(records.encode()).hexdigest(), source_files=len(files),
    dependency_versions={p.key: p.version for p in dependencies},
    dependency_locations={p.key: str(Path(p.location).resolve()) for p in dependencies})))
"""

# Both discovery helpers run under the ThorAxe lock and activated environment.
function _thoraxe_python_identity()
    paths = [Sys.which(name) for name in ("thoraxe", "transcript_query")]
    any(isnothing, paths) && error("Missing ThorAxe CLI commands.")
    length(unique(dirname.(realpath.(paths)))) == 1 ||
        error("ThorAxe CLI commands have different import directories.")
    shebangs = readline.(paths)
    all(startswith(line, "#!/") && !occursin(r"\s", line)
    for line in shebangs) || error("Unsupported ThorAxe CLI interpreter.")
    # The invocation path selects a virtual environment even when Python is a symlink.
    interpreters = [line[3:end] for line in shebangs]
    length(unique(interpreters)) == 1 ||
        error("ThorAxe CLI commands use different interpreters.")
    command = `$(first(interpreters)) -B -c $_THORAXE_IDENTITY_PROBE $(first(paths))`
    metadata = JSON.parse(read(pipeline(command; stderr = devnull), String))
    python = _executable_identity(metadata["python_executable"])
    python === nothing && error("Unknown Python executable identity.")
    return (; python, metadata,
        cli_sha256 = Dict(name => _file_sha256(path)
        for (name, path) in zip(("thoraxe", "transcript_query"), paths)))
end

function _thoraxe_aligner_identity()
    # The wrapper may replace mismatching bytes; inspect without calling get_aligner.
    isassigned(ThorAxe.ALIGNER_PATH) || error("Unknown ProGraphMSA location.")
    aligner_path = joinpath(ThorAxe.ALIGNER_PATH[], ThorAxe.ALIGNER_NAME)
    expected = ThorAxe.ALIGNER_SHA256[ThorAxe._aligner_filename()]
    _file_sha256(aligner_path) == expected ||
        error("ProGraphMSA does not match the wrapper's expected build.")
    if occursin(r"\s", aligner_path)
        isassigned(ThorAxe.ALIGNER_SHIM_PATH) ||
            error("ProGraphMSA launcher has not been created.")
        aligner_path = ThorAxe.ALIGNER_SHIM_PATH[]
    end
    aligner = _executable_identity(aligner_path)
    aligner === nothing && error("Unknown ProGraphMSA executable identity.")
    return aligner
end

function _inspect_thoraxe_identity(inspect::Function)
    try
        return inspect()
    catch err
        err isa InterruptException && rethrow()
        reason = err isa ProcessFailedException ? "Python identity probe failed." :
                 sprint(showerror, err)
        @warn "ThorAxe tool identity unavailable; affected cache reuse is disabled." reason
        return nothing
    end
end

function _thoraxe_tool_identity(; aligner::Bool = true)
    return _inspect_thoraxe_identity() do
        _with_thoraxe_lock() do
            CondaPkg.withenv() do
                python = _thoraxe_python_identity()
                # Missing ProGraphMSA must not invalidate a known transcript-query tool.
                alignment = aligner ?
                            _inspect_thoraxe_identity(_thoraxe_aligner_identity) : nothing
                return (; python, aligner = alignment)
            end
        end
    end
end

function _thoraxe_implementation_identity()
    return _implementation_identity((ThorAxe, CondaPkg, HHsuite_jll))
end

function _transcript_query_tool_identity()
    tools = _thoraxe_tool_identity(; aligner = false)
    return tools === nothing ? nothing : tools.python
end

function _thoraxe_stage_provenance(sample_count::Integer)
    tools = _thoraxe_tool_identity()
    return (; implementation = _thoraxe_implementation_identity(),
        python = tools === nothing ? nothing : tools.python,
        aligner = tools === nothing ? nothing : tools.aligner,
        hhsuite = sample_count == 0 ? nothing :
                  _jll_tool_identity(HHsuite_jll, (:hhmake, :hhalign)))
end

function _known_thoraxe_provenance(provenance, sample_count::Integer)
    provenance === nothing && return false
    all(name -> get(provenance, name, nothing) !== nothing,
        (:implementation, :python, :aligner)) || return false
    return sample_count == 0 || get(provenance, :hhsuite, nothing) !== nothing
end
