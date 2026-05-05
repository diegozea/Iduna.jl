@testitem "explicit imports lint" begin
    using Iduna
    import ExplicitImports

    ExplicitImports.test_explicit_imports(Iduna; all_qualified_accesses_are_public = false)
end
