# Repository Guidelines

## Structure

Iduna is a Julia package. Source starts at `src/Iduna.jl`; tests live in `test/runtests.jl`; docs use Documenter via `docs/make.jl` and `docs/src/index.md`; CI is in `.github/workflows/`.

## Commands

- `julia --project=. -e 'using Pkg; Pkg.test()'`: run tests.
- `julia --project=docs docs/make.jl`: build docs.

## Style

Format with JuliaFormatter. Use `CamelCase` for types/modules, `snake_case` for functions/variables, and `!` for mutating functions. Keep exports deliberate and public APIs documented.

## Tests

Group related tests into `test/*.jl` and include them from `test/runtests.jl`. Keep Aqua checks passing.

## Agent Notes

Do not use subagents by default. Use Julia MCP for short evaluation and formatting when available. Use `julia -e` for long tasks such as full tests because MCP sessions can time out. Preserve unrelated local changes.
