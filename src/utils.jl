

using ACE.Random: rand_radial, rand_sphere

function ZμRnYlm_1pbasis(; init = true, species = nothing, maxdeg = Inf, Deg = ACE.NaiveTotalDegree(), kwargs...)
   RnYlm = ACE.Utils.RnYlm_1pbasis(; init=false, D = Deg, kwargs...)
   Zμ = Species1PBasis(species)
   B1p = Zμ * RnYlm
   if init 
      ACE.init1pspec!(B1p; maxdeg = maxdeg, Deg = Deg)
   end
   return B1p
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
   X0 = AtomState( rand(Zμ), 0 * rand_sphere() )
   return AtomicEnvironment(X0, Xs)
end


