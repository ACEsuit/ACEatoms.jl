


# --------------------- Atom related States

import Base: promote_rule, *, +, real

const AtomState{T} = ACE.State{(:mu, :mu0, :rr), Tuple{AtomicNumber, AtomicNumber, SVector{3, T}}} 

+(X::TX, u::SVector{3}) where {TX <: AtomState} = TX( (rr = X.rr + u, mu = X.mu) )
+(u::SVector{3}, X::TX) where {TX <: ACE.DState{(:rr,)}} = u + X.rr

real(X::AtomState{T}) where {T} = 
      AtomState{real(T)}( (rr = real.(X.rr), mu = X.mu, mu0 = X.mu0) )
