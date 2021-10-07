import ACE: OneParticleBasis, AbstractState, Scal1pBasis, VectorPool
import ACE.OrthPolys: TransformedPolys, transformed_jacobi
import ACE.Transforms: IdTransform 
import NamedTupleTools
using NamedTupleTools: namedtuple

export Pop1PBasis

@doc raw"""
`struct Pop1PBasis <: OneParticleBasis`

One-particle basis similar to Scal1pBasis, but it has always length 1 such that P(x) = x
"""
mutable struct Pop1PBasis{TransformedPolys} <: OneParticleBasis{TransformedPolys}
  P::TransformedPolys
end

symbols(basis::Pop1PBasis) = [:μ]

indexrange(basis::Pop1PBasis) = Dict( :μ => length(basis) )

#_getidx(b, basis::Pop1PBasis) = b[_idxsym(basis) ]

isadmissible(b, basis::Pop1PBasis) = (1 <= _getidx(b, basis) <= length(basis))

function Pop1PBasis(num_basis::Int)
  P = transformed_jacobi(num_basis, IdTransform(), 0.0, 0.0)
  return Pop1PBasis(P)
end

Base.length(basis::Pop1PBasis) = length(basis.P)

_val(X::AbstractState, basis::Pop1PBasis) = X.population

function ACE.evaluate!(B, basis::Pop1PBasis, x::Number)
  return evaluate!(B, basis.P, x)
end

function ACE.evaluate!(B, basis::Pop1PBasis, X::AbstractState)
  return evaluate!(B, basis.P, _val(X, basis)) 
end

degree(b, basis::Pop1PBasis, args...) = degree(basis.P)

write_dict(V::Pop1PBasis) = 
      Dict( "__id__" => "Pop1PBasis",
      "Polys" => write_dict(V.P))

function read_dict(::Val{:Pop1PBasis}, D::Dict) 
  P = read_dict(D["Polys"])
   return Pop1PBasis(P)
end
