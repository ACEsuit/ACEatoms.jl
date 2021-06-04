


# --------------------- Atom related States

import Base: promote_rule, *, +, real

const AtomState{T} = ACE.State{(:mu, :rr), Tuple{AtomicNumber, SVector{3, T}}} 


*(A::AbstractMatrix, X::AtomState{T}) where {T} = DAtomState{T}( (rr = A * X.rr,) )

+(X::TX, u::SVector{3}) where {TX <: AtomState} = TX( (rr = X.rr + u, mu = X.mu) )
+(u::SVector{3}, X::TX) where {TX <: ACE.DState{(:rr,)}} = u + X.rr

real(X::AtomState{T}) where {T} = 
      AtomState{real(T)}( (rr = real.(X.rr), mu = X.mu) )



struct AtomicEnvironment{STT} <: AbstractConfiguration
   X0::STT   # state of center atom
   Xs::Vector{STT}  # states of neighbouring atoms
end
