

using ACE.Random: rand_radial, rand_sphere

function ZμRnYlm_1pbasis(; species = nothing, kwargs...)
   RnYlm = ACE.Utils.RnYlm_1pbasis(; kwargs...)
   Zμ = Species1PBasis(species)
   return Zμ * RnYlm
end

_rand_atstate(Zμ, Rn) = 
      AtomState( rand(Zμ), rand_radial(Rn) * rand_sphere())


function rand_environment(B1p, Nat::Integer)
   @assert Set((:μ, :n, :l, :m)) == Set(ACE.symbols(B1p))
   Rn = B1p.bases[2]
   @assert Rn isa ACE.Rn1pBasis
   Zμ = B1p.bases[1] 
   @assert Zμ isa Species1PBasis
   
   Xs = [ _rand_atstate(Zμ, Rn) for _ = 1:Nat ]
   X0 = _rand_atstate(Zμ, Rn)
   return AtomicEnvironment(X0, Xs)
end