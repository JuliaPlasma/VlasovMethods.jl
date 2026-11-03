using ExplicitImports
using GeometricMachineLearning
using Test
using VlasovMethods

# The extension loads with `GeometricMachineLearning`, so its module exists from here on.
const EXT = Base.get_extension(VlasovMethods, :MultiModalAutoEncoderExt)

@testset "ExplicitImports" begin
    @test EXT isa Module

    # The package lists no explicit import that it does not use, and takes each explicit import
    # and each qualified access from the module that owns the name.
    @test check_no_stale_explicit_imports(VlasovMethods) === nothing
    @test check_all_explicit_imports_via_owners(VlasovMethods) === nothing
    @test check_all_qualified_accesses_via_owners(VlasovMethods) === nothing

    # The multi-modal extension imports each name explicitly, from the module that owns it and
    # declares it public, and imports nothing it does not use.
    @test check_no_implicit_imports(EXT) === nothing
    @test check_no_stale_explicit_imports(EXT) === nothing
    @test check_all_explicit_imports_via_owners(EXT) === nothing
    @test check_all_qualified_accesses_via_owners(EXT) === nothing
    @test check_all_explicit_imports_are_public(EXT) === nothing
    @test check_all_qualified_accesses_are_public(EXT) === nothing
end
