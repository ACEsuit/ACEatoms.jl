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
mutable struct Pop1PBasis <: OneParticleBasis{Float64}
  # P::TransformedPolys
end

symbols(basis::Pop1PBasis) = [:P]

indexrange(basis::Pop1PBasis) = Dict( :P => 1:1 ) # Dict( :P => 1:length(basis) )

isadmissible(b, basis::Pop1PBasis) = (b.P == 1) # (1 <= b.P <= length(basis))

get_index(basis::Pop1PBasis, b) = 1

#function Pop1PBasis() #(num_basis::Int)
  #P = transformed_jacobi(num_basis, IdTransform(), 120.0, 0.0)
#  return Pop1PBasis(P)
#end

Base.length(basis::Pop1PBasis) = 1  # length(basis.P)

_val(X::AbstractState, basis::Pop1PBasis) = X.population

# Dict(:C => (4.0, 7.0), :H => (0.0, 2.0))

# (x-(a+b)/2)/(b-a)

function ACE.evaluate!(B, basis::Pop1PBasis, x::Number)
  B[1] = x 
  return B
  #return evaluate!(B, basis.P, x)
end

function ACE.evaluate!(B, basis::Pop1PBasis, X::AbstractState)
  return evaluate!(B, basis.P, _val(X, basis)) 
end

degree(b, basis::Pop1PBasis, args...) = 0 #degree(basis.P)

write_dict(V::Pop1PBasis) = 
      Dict( "__id__" => "ACE_Pop1PBasis")

function read_dict(::Val{:ACE_Pop1PBasis}, D::Dict) 
   return Pop1PBasis()
end
