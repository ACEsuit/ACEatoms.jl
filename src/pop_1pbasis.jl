import ACE: OneParticleBasis, AbstractState, Scal1pBasis
import ACE.OrthPolys: TransformedPolys
import NamedTupleTools
using NamedTupleTools: namedtuple


@doc raw"""
`struct Pop1pBasis <: OneParticleBasis`

One-particle basis similar to Scal1pBasis, but it has always length 1 such that P(x) = x
"""
mutable struct Pop1pBasis{VSYM, VIDX, ISYM, T, TT, TJ} <: Scal1pBasis{T}
   P::TransformedPolys{T, TT, TJ}
   B_pool::VectorPool{T}
   dB_pool::VectorPool{T}   
end


pop1pbasis(varsym::Symbol, idxsym::Symbol, args...; varidx = 1, kwargs...) = 
            Scal1pBasis(varsym, varidx, idxsym,  
                  ACE.OrthPolys.transformed_jacobi(args...; kwargs...)[1])

Pop1pBasis(varsym::Symbol, idxsym::Symbol, P::TransformedPolys{T, TT, TJ}
            ) where {T, TT, TJ} = Scal1pBasis(varsym, 1, idxsym, P[1])

Pop1pBasis(varsym::Symbol, varidx::Integer, idxsym::Symbol, P::TransformedPolys{T, TT, TJ}
            ) where {T, TT, TJ} = 
      Scal1pBasis{varsym, Int(varidx), idxsym, T, TT, TJ}(P[1])

Pop1pBasis{VSYM, VIDX, ISYM, T, TT, TJ}(P::TransformedPolys{T, TT, TJ}
            ) where {VSYM, VIDX, ISYM, T, TT, TJ} = 
      Scal1pBasis{VSYM, VIDX, ISYM, T, TT, TJ}(P[1], VectorPool{T}(), VectorPool{T}())

_varsym(basis::Pop1pBasis{VSYM}) where {VSYM} = VSYM
_varidx(basis::Pop1pBasis{VSYM, VIDX}) where {VSYM, VIDX} = VIDX
_idxsym(basis::Pop1pBasis{VSYM, VIDX, ISYM}) where {VSYM, VIDX, ISYM} = ISYM

_val(X::AbstractState, basis::Pop1pBasis) = 
      getproperty(X, _varsym(basis))[_varidx(basis)]

_val(x::Number, basis::Pop1pBasis) = x

# ---------------------- Implementation of Pop1pBasis


Base.length(basis::Pop1pBasis) = length(basis.P)

get_spec(basis::Pop1pBasis) =
      [  NamedTuple{(_idxsym(basis),)}(n) for n = 1:length(basis) ]

==(P1::Pop1pBasis, P2::Pop1pBasis) = (P1.P == P2.P)

write_dict(basis::Pop1pBasis{T}) where {T} = Dict(
      "__id__" => "ACE_Pop1pBasis",
          "P" => write_dict(basis.P) , 
          "varsym" => string(_varsym(basis)), 
          "varidx" => _varidx(basis), 
          "idxsym" => string(_idxsym(basis)) )

read_dict(::Val{:ACE_Pop1pBasis}, D::Dict) =   
      Pop1pBasis(Symbol(D["varsym"]), Int(D["varidx"]), Symbol(D["idxsym"]), 
                  read_dict(D["P"]))

valtype(basis::Pop1pBasis) = valtype(basis.P)

valtype(basis::Pop1pBasis, cfg::AbstractConfiguration) = 
      valtype(basis, zero(eltype(cfg)))

valtype(basis::Pop1pBasis, X::AbstractState) = 
      valtype(basis.P, _val(X, basis))

gradtype(basis::Pop1pBasis, X::AbstractState) = 
      dstate_type(valtype(basis, X), X)

symbols(basis::Pop1pBasis) = [ _idxsym(basis) ]

indexrange(basis::Pop1pBasis) = NamedTuple{(_idxsym(basis), )}((1:length(basis),))

_getidx(b, basis::Pop1pBasis) = b[_idxsym(basis) ]

isadmissible(b, basis::Pop1pBasis) = (_getidx(b, basis) == 1)

degree(b, basis::Pop1pBasis) = _getidx(b, basis) - 1

get_index(basis::Pop1pBasis, b) = _getidx(b, basis)

# TODO: need better structure to support this ... 
rand_radial(basis::Pop1pBasis) = rand_radial(basis.P)

# ---------------------------  Evaluation code
#

evaluate!(B, basis::Pop1pBasis, x::Number) =
      evaluate!(B, basis.P, x)

evaluate!(B, basis::Pop1pBasis, X::AbstractState) =
      evaluate!(B, basis.P, _val(X, basis))

"""
returns an `SVector{N}` of the form `x * e_I` where `e_I` is the Ith canonical basis vector.
"""
@generated function __e(::SVector{N}, ::Val{I}, x::T) where {N, I, T}
   code = "SA["
   for i = 1:N 
      if i == I
         code *= "x,"
      else 
         code *= "0,"
      end
   end
   code *= "]"
   quote 
      $( Meta.parse(code) )
   end
end

__e(::Number, ::Any, x) = x


function _pop1pbasis_grad(TDX::Type, basis::Pop1pBasis, gval)
   gval_tdx = __e( getproperty(zero(TDX), _varsym(basis)), 
                   Val(_varidx(basis)), 
                   gval )
   return TDX( NamedTuple{(_varsym(basis),)}((gval_tdx,)) )
end

function evaluate_d!(dB, basis::Pop1pBasis, X::AbstractState)
   TDX = eltype(dB)
   x = _val(X, basis)
   dP = acquire_dB!(basis.P, x)
   evaluate_d!(dP, basis.P, x)
   for n = 1:length(basis)
      dB[n] = _pop1pbasis_grad(TDX, basis, dP[n])
   end
   release_dB!(basis.P, dP)
   return dB
end

function evaluate_ed!(B, dB, basis::Pop1pBasis, X::AbstractState)
   TDX = eltype(dB)
   x = _val(X, basis)
   dP = acquire_dB!(basis.P, x)
   evaluate!(B, basis.P, x)
   evaluate_d!(dP, basis.P, x)
   for n = 1:length(basis)
      dB[n] = _pop1pbasis_grad(TDX, basis, dP[n])
   end
   release_dB!(basis.P, dP)
   return B, dB
end


# this one we probably only need for training so can relax the efficiency a bit 
function evaluate_dd(basis::Pop1pBasis, X::AbstractState) 
   ddP = ForwardDiff.derivative(x -> evaluate_d(basis, _val(X, basis)))
   TDX = gradtype(basis, X)
   return _pop1pbasis_grad.(Ref(TDX), Ref(basis), ddP_n)
end



# -------------- AD codes 

import ChainRules: rrule, ZeroTangent, NoTangent

function _rrule_evaluate(basis::Pop1pBasis, X::AbstractState, 
                         w::AbstractVector{<: Number})
   @assert _varidx(basis) == 1
   x = _val(X, basis)
   a = _rrule_evaluate(basis.P, x, real.(w))
   TDX = ACE.dstate_type(a, X)
   return TDX( NamedTuple{(_varsym(basis),)}( (a,) ) )
end

rrule(::typeof(evaluate), basis::Pop1pBasis, X::AbstractState) = 
                  evaluate(basis, X), 
                  w -> (NoTangent(), NoTangent(), _rrule_evaluate(basis, X, w))

             
                  
function _rrule_evaluate_d(basis::Pop1pBasis, X::AbstractState, 
                           w::AbstractVector)
   @assert _varidx(basis) == 1
   x = _val(X, basis)
   w1 = [ _val(w, basis) for w in w ]
   a = _rrule_evaluate_d(basis.P, x, w1)
   TDX = ACE.dstate_type(a, X)
   return TDX( NamedTuple{(_varsym(basis),)}( (a,) ) )
end

function rrule(::typeof(evaluate_d), basis::Pop1pBasis, X::AbstractState)
   @assert _varidx(basis) == 1
   x = _val(X, basis)
   dB_ = evaluate_d(basis.P, x)
   TDX = dstate_type(valtype(basis, X), X)
   dB = [ TDX( NamedTuple{(_varsym(basis),)}( (dx,) ) )  for dx in dB_ ]
   return dB, 
          w -> (NoTangent(), NoTangent(), _rrule_evaluate_d(basis, X, w))
end
