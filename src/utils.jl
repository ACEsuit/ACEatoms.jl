
using ACE.Random: rand_radial, rand_sphere
using ACE: ACEConfig, State 

"""
`ZμRnYlm_1pbasis` : utility function to quickly generate a 
`Zμ * Rn * Ylmn` 1-particle basis.

todo - finish the docs.
"""
function ZμRnYlm_1pbasis(; init = true, species = nothing, maxdeg = nothing, 
                           maxL = maxdeg, 
                           Bsel = ACE.SimpleSparseBasis(1, maxdeg), 
                           kwargs...)
   RnYlm = ACE.Utils.RnYlm_1pbasis(; maxdeg=maxdeg, maxL=maxL, Bsel = Bsel, kwargs...)
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


function rand_ACEConfig(B1p, Nat::Integer)
   @assert Set((:μ, :n, :l, :m)) == Set(ACE.symbols(B1p))
   Rn = B1p["Rn"]
   Zμ = B1p.bases[1] 
   @assert Zμ isa Species1PBasis
   
   mu0 = rand(Zμ)
   Xs = [ _rand_atstate(mu0, Zμ, Rn) for _ = 1:Nat ]
   return ACEConfig(Xs)
end


"""
`PopZμRnYlm_1pbasis` : utility function to quickly generate a 
`Pop * Zμ * Rn * Ylmn` 1-particle basis.
"""
function PopZμRnYlm_1pbasis(; init = true, species = nothing, maxdeg = nothing, 
                           maxL = maxdeg, 
                           Bsel = ACE.SimpleSparseBasis(1, maxdeg), 
                           kwargs...)
   RnYlm = ACE.Utils.RnYlm_1pbasis(; maxdeg=maxdeg, maxL=maxL, Bsel = Bsel, kwargs...)
   Zμ = Species1PBasis(species)
   Pop = Pop1PBasis()
   B1p = ACE.Product1pBasis( (Pop, Zμ, RnYlm.bases...) )
   if init 
      ACE.init1pspec!(B1p, Bsel)
   end
   return B1p
end

_rand_atstate(mu0, Zμ, Rn, Pop) = 
     PopAtomState{Float64}(mu = rand(Zμ), mu0 = mu0, 
                  rr = rand_radial(Rn) * rand_sphere(),
                  pop = rand() )

function rand_ACEConfig_pop(B1p, Nat::Integer)
   @assert Set((:μ, :n, :l, :m, :P)) == Set(ACE.symbols(B1p))
   Rn = B1p["Rn"]
   Zμ = B1p.bases[2] 
   @assert Zμ isa Species1PBasis
   Pop = B1p["Pop"]
   
   mu0 = rand(Zμ)
   Xs = [ _rand_atstate(mu0, Zμ, Rn, Pop) for _ = 1:Nat ]
   return ACEConfig(Xs)
end

function pair_basis(; species = :X,
   # transform parameters
   r0 = 2.5,
   trans = PolyTransform(2, r0),
   # degree parameters
   maxdeg = 8,
   # radial basis parameters
   rcut = 5.0,
   rin = 0.5 * r0,
   pcut = 2,
   pin = 0,
   rbasis = transformed_jacobi(maxdeg, trans, rcut, rin; pcut=pcut, pin=pin))
   return PolyPairBasis(rbasis, species)
end
