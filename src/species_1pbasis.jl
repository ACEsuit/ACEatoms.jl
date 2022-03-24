
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

"""
`Species1PBasis` : todo write docs 
"""
struct Species1PBasis{NZ} <: AbstractSpecies1PBasis{NZ}
   zlist::SZList{NZ}
end

zlist(basis::AbstractSpecies1PBasis) = basis.zlist

Species1PBasis(species) = Species1PBasis(ZList(species, static=true))

Base.length(basis::AbstractSpecies1PBasis) = length(basis.zlist)

function ACE.evaluate!(A, basis::Species1PBasis, X::AbstractState)
   fill!(A, false)
   A[z2i(basis.zlist, X.mu)] = true
   return A
end

ACE.valtype(::AbstractSpecies1PBasis, args...) = Bool

symbols(::Species1PBasis) = [:μ]

indexrange(basis::Species1PBasis) = Dict( :μ => Int.(basis.zlist.list) )

isadmissible(b, basis::Species1PBasis) = (AtomicNumber(b.μ) in basis.zlist.list)

degree(b, basis::AbstractSpecies1PBasis, args...) = 0

get_index(basis::Species1PBasis, b) = z2i(basis.zlist, AtomicNumber(b.μ))

ACE.get_spec(basis::Species1PBasis, i::Integer) = (μ = basis.zlist.list[i],)
ACE.get_spec(basis::Species1PBasis) = ACE.get_spec.(Ref(basis), 1:length(basis))

ACE._scaling_absvalue(x::AtomicNumber) = 0


Base.rand(basis::AbstractSpecies1PBasis) = rand(basis.zlist.list)

write_dict(V::Species1PBasis) = 
      Dict( "__id__" => "Species1PBasis", 
            "zlist" => write_dict(V.zlist))

function read_dict(::Val{:Species1PBasis}, D::Dict) 
   zlist = read_dict(D["zlist"])
   return Species1PBasis(zlist)
end
