
import ACE

import ACE: OneParticleBasis,
            symbols,
            indexrange,
            isadmissible,
            degree,
            get_index

import ACEbase: Discrete1pBasis

import JuLIP.Potentials: ZList, SZList, z2i, i2z


export Species1PBasis


"""
and auxiliary type to write a joint codebase for
   `Species1PBasis` and `Species1PBasisCtr`
"""
abstract type AbstractSpecies1PBasis{NZ} <: Discrete1pBasis{Bool} end

struct Species1PBasis{NZ} <: AbstractSpecies1PBasis{NZ}
   zlist::SZList{NZ}
end

zlist(basis::AbstractSpecies1PBasis) = basis.zlist

Species1PBasis(species) = Species1PBasis(ZList(species, static=true))

Base.length(basis::AbstractSpecies1PBasis) = length(basis.zlist)^2

function ACE.evaluate!(A, basis::Species1PBasis, X::AbstractState)
   fill!(A, false)
   B[z2i(basis.zlist, X.mu)] = true
   return B
end

ACE.valtype(::AbstractSpecies1PBasis, X) = Bool

symbols(::Species1PBasis) = [:μ]

indexrange(basis::Species1PBasis) = Dict( :μ => Int.(basis.zlist.list) )

isadmissible(b, basis::Species1PBasis) = (AtomicNumber(b.μ) in basis.zlist.list)

degree(b, basis::AbstractSpecies1PBasis, args...) = 0

get_index(basis::Species1PBasis, b) = z2i(basis.zlist, AtomicNumber(b.μ))

Base.rand(basis::AbstractSpecies1PBasis) = rand(basis.zlist.list)