
import ACE

import ACE: OneParticleBasis,
            symbols,
            indexrange,
            isadmissible,
            degree,
            get_index

import ACEbase: Discrete1pBasis

import JuLIP.Potentials: ZList, SZList, z2i, i2z


export Species1PBasis, Species1PBasisCtr


"""
and auxiliary type to write a joint codebase for
   `Species1PBasis` and `Species1PBasisCtr`
"""
abstract type AbstractSpecies1PBasis{NZ} <: Discrete1pBasis{Bool} end

struct Species1PBasis{NZ} <: AbstractSpecies1PBasis{NZ}
   zlist::SZList{NZ}
end

struct Species1PBasisCtr{NZ} <: AbstractSpecies1PBasis{NZ}
   zlist::SZList{NZ}
end


zlist(basis::AbstractSpecies1PBasis) = basis.zlist

Species1PBasis(species) = Species1PBasis(ZList(species, static=true))
Species1PBasisCtr(species) = Species1PBasisCtr(ZList(species, static=true))

Base.length(basis::AbstractSpecies1PBasis) = length(basis.zlist)^2

# ACE.evaluate!(B, tmp, basis::Species1PBasisCtr,
#               cfg::AtomicEnvironment) = evaluate!(B, tmp, basis, cfg.X0)
# ACE.evaluate!(B, tmp, basis::Species1PBasis,
#               cfg::AtomicEnvironment) = evaluate!(B, tmp, basis, Xj)

function ACE.evaluate!(B, tmp, basis::Species1PBasis, X::AbstractState)
   fill!(B, 0)
   B[z2i(basis.zlist, X.mu)] = 1
   return B
end

# technically, this is the interface, but the Species1PBasis should never be
# used outside the context of a product basis, and there it will always 
# be called via ACE.evaluate!
# function ACE.add_into_A!(A, tmp, basis::Species1PBasis{NZ}, X::AbstractState) where {NZ}
#    A[z2i(basis.zlist, X.mu)] += 1
#    return nothing
# end

ACE.fltype(::AbstractSpecies1PBasis) = Bool

ACE.gradtype(basis::Species1PBasis, X::AbstractState) = 
         ACE.dstate_type(fltype(basis), X)

symbols(::Species1PBasis) = [:μ]
symbols(::Species1PBasisCtr) = [:μ0]

indexrange(basis::Species1PBasis) = Dict( :μ => Int.(basis.zlist.list) )
indexrange(basis::Species1PBasisCtr) = Dict( :μ0 => Int.(basis.zlist.list) )

isadmissible(b, basis::Species1PBasis) = (AtomicNumber(b.μ) in basis.zlist.list)
isadmissible(b, basis::Species1PBasisCtr) = (AtomicNumber(b.μ0) in basis.zlist.list)

degree(b, basis::AbstractSpecies1PBasis, args...) = 0

get_index(basis::Species1PBasis, b) = z2i(basis.zlist, AtomicNumber(b.μ))
get_index(basis::Species1PBasisCtr, b) = z2i(basis.zlist, AtomicNumber(b.μ0))

Base.rand(basis::AbstractSpecies1PBasis) = rand(basis.zlist.list)