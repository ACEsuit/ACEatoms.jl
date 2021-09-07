

using ACE.Random: rand_radial, rand_sphere
using ACE: ACEConfig

function ZμRnYlm_1pbasis(; init = true, species = nothing, maxdeg = nothing, 
                           maxL = maxdeg, 
                           Bsel = ACE.SimpleSparseBasis(1, maxdeg), kwargs...)
   RnYlm = ACE.Utils.RnYlm_1pbasis(; maxdeg=maxdeg, maxL=maxL, Bsel = Bsel)
   Zμ = Species1PBasis(species)
   B1p = Zμ * RnYlm
   if init 
      ACE.init1pspec!(B1p, Bsel)
   end
   return B1p
end

_rand_atstate(mu0, Zμ, Rn) = 
      AtomState{Float64}(mu = rand(Zμ), 
                         mu0 = mu0, 
                         rr = rand_radial(Rn) * rand_sphere() )


function rand_environment(B1p, Nat::Integer)
   @assert Set((:μ, :n, :l, :m)) == Set(ACE.symbols(B1p))
   Rn = B1p.bases[2]
   @assert Rn isa ACE.Rn1pBasis
   Zμ = B1p.bases[1] 
   @assert Zμ isa Species1PBasis
   
   mu0 = rand(Zμ)
   Xs = [ _rand_atstate(mu0, Zμ, Rn) for _ = 1:Nat ]
   return ACEConfig(Xs)
end


