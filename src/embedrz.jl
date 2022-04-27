

import ACE, ACEbase

import ACE: SList, i2val, val2i, get_spec, symbols, indexrange, get_index, degree 

using ACE: @λ, chain, GetNorm, B1pComponent

# -------------------------


struct EmbedRZ{RSYM, SSYM, TPR, T}
   Pr::TPR
   embeddings::Dict{Any, Matrix{T}}
end

function EmbedRZ(rsym::Symbol, ssym::Symbol, Pr::TPR, embeddings::Dict
                )  where {TPR}
   T = eltype(first(embeddings)[2])
   embeddings = Dict{Any, Matrix{T}}(embeddings)
   return EmbedRZ{rsym, ssym, TPR, T}(Pr, embeddings)
end


_rsym(basis::EmbedRZ{RSYM, SSYM}) where {RSYM, SSYM} = RSYM
_ssym(basis::EmbedRZ{RSYM, SSYM}) where {RSYM, SSYM} = SSYM

_rr(basis::EmbedRZ, X) = getproperty(X, _rsym(basis))
_μ(basis::EmbedRZ, X) = getproperty(X, _ssym(basis))


function embed_nl_spec1(maxdeg::Integer, wl = 1.5)
   spec = [] 
   for n = 1:maxdeg
      for l = 0:floor(Int, (maxdeg - n)/wl)
         push!(spec, (n = n, l = l))
      end
   end 
   return identity.(spec)
end

embed_nl_spec_max(maxn::Integer, maxl::Integer) = 
      [ (n = n, l = l) for n = 1:maxn, l = 0:maxl ]

function embedrz(Pr, 
                 spec::AbstractVector{<: NamedTuple}, 
                 embeddings::Dict; 
                 rrsym = :rr, ssym = :mu, label = "Rnl")
   getval = GetNorm{rrsym}() 
   degrees = zeros(Int, length(spec))
   Rnl = EmbedRZ(rrsym, ssym, Pr, embeddings)
   return B1pComponent(Rnl, getval, spec, degrees, label)
end


function evaluate(basis::EmbedRZ, X)
   P = evaluate(basis.Pr, _rr(basis, X))
   uP = basis.embeddings[_μ(basis, X)] * P 
   ACE.release!(P)
   return uP 
end


