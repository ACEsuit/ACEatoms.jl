using ACEatoms
using Test

@testset "ACEatoms.jl" begin

   # ----------------------
   #   pair potentials
   include("pair/test_pair_basis.jl")
   include("pair/test_pair_pot.jl")
   include("pair/test_repulsion.jl")

   # ----------------------
   #   species 1p basis 
   @testset "ACEatoms 1p Basis" include("test_1pbasis.jl") 
   #   basic site energy calculator
   include("test_siteenergy.jl")

   # ---------------------- 
   #  special physics 
   include("test_electro.jl")
end
