


# --------------------- Atom related States


"same as EuclideanVectorState"
const PositionState{T} = EuclideanVectorState{T}

@doc raw"""
`struct SpeciesState` : a $\mathbb{Z}$ value, which is invariant under the
rotation group. It defines an atomic species.
"""
struct SpeciesState <: AbstractDiscreteState
   mu::AtomicNumber
   label::String
   SpeciesState(z_or_sym, label::String = "") =
         new(AtomicNumber(z_or_sym), label)
end

SpeciesState(label::String = "") = SpeciesState(AtomicNumber(0), label)

Base.show(io::IO, s::SpeciesState) =
      print(io, "$(s.label)[$(chemical_symbol(s.mu))]")



"""
`struct AtomState` : basic implementation of the state of an Atom
consistent with original ACE.
"""
struct AtomState{T} <: AbstractState
   mu::AtomicNumber
   rr::SVector{3, T}
   # add other features here? 
   # - charge 
   # - dipole 
   # ...
end

AtomState(mu, rr::AbstractVector{T}) where {T} =
   AtomState(AtomicNumber(mu), SVector{3, T}(rr...))
AtomState(T::Type) = AtomState(0, zero(SVector{3, T}))
AtomState(mu = 0, T::Type = Float64) = AtomState(mu, zero(SVector{3, T}))

Base.show(io::IO, X::AtomState) =
   print(io, ( SpeciesState(X.mu, "μ"),
               PositionState(X.rr, "𝒓") ))


struct AtomicEnvironment{STT} <: AbstractConfiguration
   X0::STT   # state of center atom
   Xs::Vector{STT}  # states of neighbouring atoms
end
