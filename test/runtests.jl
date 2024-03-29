using ACEatoms
using Test
using JuLIP 
JuLIP._usethreads[] = false 

##

@testset "ACEatoms.jl" begin

   # ----------------------
   #   pair potentials
   @testset "PolyPairBasis" begin include("pair/test_pair_basis.jl") end 
   @testset "PolyPairPot" begin include("pair/test_pair_pot.jl") end 
   @testset "RepulsiveCore" begin include("pair/test_repulsion.jl") end 

   # ----------------------
   #   species 1p basis 
   @testset "ACEatoms 1p Basis" begin include("test_1pbasis.jl")  end 
   #   basic site energy calculator
   @testset "ACESitePotential" begin  include("test_siteenergy.jl") end 

   # ---------------------- 
   #  special physics 
   @testset "Electrostatics" begin include("test_electro.jl") end 
   @testset "Dipoles" begin include("test_dip.jl") end 


   @testset "z Embedding" begin include("test_embed1p.jl") end 
   @testset "(r,z) Embedding" begin include("test_embedrz.jl") end 

end
