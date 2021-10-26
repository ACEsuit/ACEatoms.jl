

using ACE.Random: rand_radial, rand_sphere
using ACE: ACEConfig, State 
using ACE.Utils: BondBasisSelector, RnYlm_1pbasis
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
   Rn = B1p.bases[2]
   @assert Rn isa ACE.Rn1pBasis
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
   B1p = Pop * Zμ * RnYlm
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
   Rn = B1p.bases[3]
   @assert Rn isa ACE.Rn1pBasis
   Zμ = B1p.bases[2] 
   @assert Zμ isa Species1PBasis
   Pop = B1p.bases[1] 
   @assert Pop isa Pop1PBasis
   
   mu0 = rand(Zμ)
   Xs = [ _rand_atstate(mu0, Zμ, Rn, Pop) for _ = 1:Nat ]
   return ACEConfig(Xs)
end


#isym=:be, bond_weight = 1.0, env_weight = 1.0)
function SymmetricBondSpecies_basis(ϕ::ACE.AbstractProperty, env::ACE.BondEnvelope, Bsel::ACE.SparseBasis; species = nothing, RnYlm = nothing, kwargs...)
   BondSelector =  BondBasisSelector(Bsel; kwargs...)
   if RnYlm === nothing
       RnYlm = RnYlm_1pbasis(;   r0 = ACE.cutoff_radialbasis(env), 
                                           rin = 0.0,
                                           trans = PolyTransform(2, ACE.cutoff_radialbasis(env)), 
                                           pcut = 2,
                                           pin = 0, 
                                           kwargs...
                                       )
   end
   Zμ = Species1PBasis(species)
   Bc = ACE.Categorical1pBasis([:bond, :env]; varsym = :be, idxsym = :be )
   B1p = Zμ * Bc * RnYlm * env
   return ACE.SymmetricBasis(ϕ, B1p, BondSelector)
end
