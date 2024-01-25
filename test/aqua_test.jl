using FunctionZeros
using Aqua: Aqua

const PkgName = FunctionZeros

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(PkgName)
end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(PkgName)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(PkgName)
end

if VERSION >= v"1.1"
    # Depending on Optim causes many ambiguity errors outside our control
    @testset "aqua test ambiguities" begin
        Aqua.test_ambiguities([PkgName, Core, Base])
    end
end

@testset "aqua piracy" begin
    Aqua.test_piracies(PkgName)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(PkgName)
end

@testset "aqua state deps" begin
    Aqua.test_stale_deps(PkgName)
end
