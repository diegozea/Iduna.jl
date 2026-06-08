# Repository Guidelines

## Structure

Iduna is a Julia package. Source starts at `src/Iduna.jl`; tests live in
`test/runtests.jl`; docs use Documenter via `docs/make.jl` and
`docs/src/index.md`; CI is in `.github/workflows/`.

## Commands

- `julia --project=. -e 'using Pkg; Pkg.test()'`: run tests.
- `julia --project=docs docs/make.jl`: build docs.

## Style

Format with JuliaFormatter. Use `CamelCase` for types/modules, `snake_case` for
functions/variables, and `!` for mutating functions. Keep exports deliberate and
public APIs documented.

Docstrings must follow the SciMLStyle documentation guidelines: use concise
Markdown, wrap lines at 92 characters, document exported functions and public
types, and prefer the standard type/function sections (`# Fields`, `# Arguments`,
`# Keywords`, `# Returns`, `# Throws`) when they help. Write docstrings and
documentation in simple plain English. Avoid computer-science jargon unless the
term is necessary, and prefer wording that is clear to readers coming from life
sciences.

## Tests

Group related tests into `test/*.jl` and include them from `test/runtests.jl`.
Keep Aqua checks passing. Every new or changed source behavior must be covered
by tests. Do not create a PR that decreases test coverage or leaves newly added
executable source lines uncovered; before handoff, run relevant tests and a
coverage check for source changes.

## Agent Notes

Do not use subagents by default. Use Julia MCP for short evaluation and
formatting when available. Use `julia -e` for long tasks such as full tests
because MCP sessions can time out. Preserve unrelated local changes.
