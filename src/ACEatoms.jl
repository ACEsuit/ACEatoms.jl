
module ACEatoms

using StaticArrays, Reexport

include("ext_imports.jl")
import JuLIP: dipole, energy, forces

include("configs.jl")

include("species_1pbasis.jl")

include("pop_1pbasis.jl")

include("pairpots/pair.jl")
@reexport using ACEatoms.PairPotentials

include("electrostatics.jl")
@reexport using ACEatoms.Electrostatics

include("utils.jl")

include("siteenergy.jl")

# DONT UNCOMMENT ANYTHING BELOW HERE FOR NOW!!!
# include("ad.jl")

# include("models/models.jl")

include("dipole.jl")

include("embedrz.jl")

end
