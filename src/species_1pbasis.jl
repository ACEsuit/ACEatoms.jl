
import ACE

import ACE: OneParticleBasis,
            symbols,
            indexrange,
            isadmissible,
            degree,
            get_index

import JuLIP.Potentials: ZList, SZList, z2i, i2z


export Species1PBasisCtr, Species1PBasisNeig


"""
and auxiliary type to write a joint codebase for
   `Species1PBasisCtr` and `Species1PBasisNeig`
"""
abstract type SpeciesBasis{NZ} <: OneParticleBasis{Bool} end

struct Species1PBasisCtr{NZ} <: SpeciesBasis{NZ}
   zlist::SZList{NZ}
end

struct Species1PBasisNeig{NZ} <: SpeciesBasis{NZ}
   zlist::SZList{NZ}
end


zlist(basis::SpeciesBasis) = basis.zlist

Species1PBasisCtr(species) = Species1PBasisCtr(ZList(species, static=true))
Species1PBasisNeig(species) = Species1PBasisNeig(ZList(species, static=true))

Base.length(basis::SpeciesBasis) = length(basis.zlist)^2

ACE.evaluate!(B, tmp, basis::Species1PBasisCtr,
              cfg::AtomicEnvironment) = evaluate!(B, tmp, basis, cfg.X0)
ACE.evaluate!(B, tmp, basis::Species1PBasisNeig,
              cfg::AtomicEnvironment) = evaluate!(B, tmp, basis, Xj)

function ACE.evaluate!(B, tmp, basis::SpeciesBasis{NZ}, X) where {NZ}
   fill!(B, 0)
   B[z2i(basis.zlist, X.mu)] = 1
   return B
end

ACE.fltype(::SpeciesBasis) = Bool

symbols(::Species1PBasisCtr) = [:μ0]
symbols(::Species1PBasisNeig) = [:μ]

indexrange(basis::Species1PBasisCtr) = Dict( :μ0 => Int.(basis.zlist.list) )
indexrange(basis::Species1PBasisNeig) = Dict( :μ => Int.(basis.zlist.list) )

isadmissible(b, basis::Species1PBasisCtr) = (AtomicNumber(b.μ0) in basis.zlist.list)
isadmissible(b, basis::Species1PBasisNeig) = (AtomicNumber(b.μ) in basis.zlist.list)

degree(b, basis::SpeciesBasis) = 0

get_index(basis::Species1PBasisCtr, b) = z2i(basis.zlist, AtomicNumber(b.μ0))
get_index(basis::Species1PBasisNeig, b) = z2i(basis.zlist, AtomicNumber(b.μ))
