
module ACEatoms

using StaticArrays

using ACE: AbstractState,
           AbstractDiscreteState,
           AbstractContinuousState,
           AbstractConfiguration,
           EuclideanVectorState,
           PolyTransform,
           transformed_jacobi

export SpeciesState, PositionState

include("julip_imports.jl")

include("configs.jl")

include("utils.jl")

include("species_1pbasis.jl")

include("pairpots/pair.jl")

end
