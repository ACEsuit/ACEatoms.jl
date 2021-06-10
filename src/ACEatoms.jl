
module ACEatoms

using StaticArrays, Reexport

include("ext_imports.jl")

using ACE: EuclideanVectorState,
           PolyTransform,
           transformed_jacobi

export SpeciesState, PositionState

import JuLIP: dipole, energy, forces

include("configs.jl")

include("species_1pbasis.jl")

include("pairpots/pair.jl")
@reexport using ACEatoms.PairPotentials

include("electrostatics.jl")
@reexport using ACEatoms.Electrostatics

include("utils.jl")

include("siteenergy.jl")

# include("ad.jl")

include("dipole.jl")

end
