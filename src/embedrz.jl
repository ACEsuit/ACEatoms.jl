

import ACE, ACEbase

import ACE: SList, i2val, val2i, get_spec, symbols, indexrange, get_index, degree, 
            evaluate, evaluate_d, evaluate_ed 

using ACE: @λ, chain, GetNorm, B1pComponent

# ---------------------------------------------------------------------------
# this is a bit of a nasty hack, still need to figure out how to do this 
# in a more systematic way

import ACE.Transforms: StaticGet, get_symbols, getval, getval_d, grad_type_dP
import ACE: rrule_evaluate, DState
using LinearAlgebra: norm 

struct GetNormMu{RSYM, SSYM} <: StaticGet end 

getval(X, ::GetNormMu{RSYM, SSYM}) where {RSYM, SSYM} = 
      (norm(getproperty(X, RSYM)), getproperty(X, SSYM))

function getval_d(X, ::GetNormMu{RSYM}) where {RSYM}
   x = getproperty(X, RSYM)
   return DState( NamedTuple{(RSYM,)}( (x/norm(x),) ) )
end 

get_symbols(::GetNormMu{RSYM, SSYM}) where {RSYM, SSYM} = (RSYM, SSYM)

function rrule_evaluate(dP, _val::GetNormMu{RSYM}, X) where {RSYM}
   rr = getproperty(X, RSYM)
   drr = rr/norm(rr)
   TDX = ACE.dstate_type(X)
   tdx = g -> TDX(DState( NamedTuple{(RSYM,)}( (g,) ) ))
   return [ tdx(drr * dP[n]) for n = 1:length(dP) ]
end


# ---------------------------------------------------------------------------


struct EmbedRZ{TPR, T}
   Pr::TPR
   embeddings::Dict{Any, Matrix{T}}
end

function EmbedRZ(Pr::TPR, embeddings::Dict)  where {TPR}
   T = eltype(first(embeddings)[2])
   embeddings = Dict{Any, Matrix{T}}(embeddings)
   return EmbedRZ(Pr, embeddings)
end


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
   getval = GetNormMu{rrsym, ssym}() 
   degrees = zeros(Int, length(spec))
   Rnl = EmbedRZ(Pr, embeddings)
   return B1pComponent(Rnl, getval, spec, degrees, label)
end


function evaluate(basis::EmbedRZ, r_μ::Tuple{<: AbstractFloat, Symbol})
   r, μ = r_μ
   P = evaluate(basis.Pr, r)
   uP = basis.embeddings[μ] * P 
   ACE.release!(P)
   return uP 
end


function evaluate_ed(basis::EmbedRZ, r_μ::Tuple{<: AbstractFloat, Symbol})
   r, μ = r_μ
   P, dP = evaluate_ed(basis.Pr, r)
   uP = basis.embeddings[μ] * P 
   udP = basis.embeddings[μ] * dP 
   ACE.release!(P)
   ACE.release!(dP)
   return uP, udP 
end

