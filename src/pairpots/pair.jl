


module PairPotentials

include("../ext_imports.jl")

using LinearAlgebra: norm, dot
using JuLIP.Potentials: ZList, SZList, @pot, @D,
                        PairPotential, SimplePairPotential
using StaticArrays: SMatrix


include("pair_basis.jl")

include("pair_pot.jl")

include("repulsion.jl")

end
