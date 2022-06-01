

import ACE, ACEbase

import ACE: SList, i2val, val2i, get_spec, symbols, indexrange, get_index, degree 


# -------------------------

@doc raw"""
`EmbedCat1pBasis` : defines an embedding of a categorical basis; 
todo - document the details... 
"""
struct EmbedCat1pBasis{VSYM, ISYM, LEN, TCAT, T} <: Discrete1pBasis{T}
   categories::SList{LEN, TCAT}
   embed::Vector{Vector{T}}
   label::String
end

_varsym(::EmbedCat1pBasis{VSYM, ISYM}) where {VSYM, ISYM} = VSYM
_isym(::EmbedCat1pBasis{VSYM, ISYM}) where {VSYM, ISYM} = ISYM

_val(X, B::EmbedCat1pBasis) = getproperty(X, _varsym(B))
_idx(b, B::EmbedCat1pBasis) = getproperty(b, _isym(B))

Base.length(B::EmbedCat1pBasis) = length( first(B.embed) )

function EmbedCat1pBasis(categories::AbstractArray, embeddings::Dict; 
              varsym::Symbol = :mu, idxsym::Symbol = :k, 
              label = "EC$(idxsym)") 
   list = SList(categories)
   lenb = length( first(embeddings)[2] )
   @assert all(lenb == length(x) for x in values(embeddings))
   embed = [ embeddings[i2val(list, i)] for i = 1:length(categories) ]
   return EmbedCat1pBasis(list, embed, varsym, idxsym, label)
end

EmbedCat1pBasis(categories::SList{LEN, TCAT}, embed::Vector{Vector{T}}, 
                varsym::Symbol, isym::Symbol, label::String) where {LEN, TCAT, T} = 
      EmbedCat1pBasis{varsym, isym, LEN, TCAT, T}(categories, embed, label)

evaluate(basis::EmbedCat1pBasis, X::AbstractState) = 
      evaluate!(Vector{_valtype(basis)}(undef, length(basis)), basis, X)


function ACE.evaluate!(A, basis::EmbedCat1pBasis, X::AbstractState)
   key = _val(X, basis)
   idx = val2i(basis.categories, key)
   copy!(A, basis.embed[idx])
   return A
end

_valtype(::EmbedCat1pBasis{VSYM, ISYM, LEN, TCAT, T}, args...
           ) where {VSYM, ISYM, LEN, TCAT, T} = T

symbols(basis::EmbedCat1pBasis) = [ _isym(basis), ]

indexrange(basis::EmbedCat1pBasis) = Dict( _isym(basis) => 1:length(basis) )

isadmissible(b, basis::EmbedCat1pBasis) = (1 <= _idx(b, basis) <= length(basis))

get_index(B::EmbedCat1pBasis, b) = _idx(b, B)

degree(b, basis::EmbedCat1pBasis, args...) = 0

Base.rand(basis::EmbedCat1pBasis) = rand(basis.list)

function get_spec(basis::EmbedCat1pBasis, i)
   return NamedTuple{(_isym(basis),)}((i,))
end

get_spec(basis::EmbedCat1pBasis) = [ get_spec(basis, i) for i = 1:length(basis) ]


write_dict(B::EmbedCat1pBasis) = 
      Dict( "__id__" => "ACE_EmbedCat1pBasis", 
            "categories" => write_dict(B.categories), 
            "T" => write_dict(valtype(B)), 
            "embed" => B.embed, 
            "VSYM" => String(_varsym(B)), 
            "ISYM" => String(_isym(B)), 
            "label" => B.label)

function read_dict(::Val{:ACE_EmbedCat1pBasis}, D::Dict)  
   T = read_dict(D["T"])
   return EmbedCat1pBasis(    
                  read_dict(D["categories"]), 
                  Vector{T}.(D["embed"]), 
                  Symbol(D["VSYM"]), 
                  Symbol(D["ISYM"]), 
                  D["label"] )
end
