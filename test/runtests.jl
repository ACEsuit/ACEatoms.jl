using ACEatoms
using Test

@testset "ACEatoms.jl" begin

    # ----------------------
    #   pair potentials
    include("pair/test_pair_basis.jl")
    include("pair/test_pair_pot.jl")
    include("pair/test_repulsion.jl")

end
