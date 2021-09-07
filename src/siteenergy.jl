

import JuLIP.Potentials: SitePotential

# TODO: nasty hack - must fix this!!!

cutoff(model::ACE.LinearACEModel) = cutoff(model.basis)

cutoff(basis::ACE.SymmetricBasis) = cutoff(basis.pibasis)

cutoff(basis::ACE.PIBasis) = cutoff(basis.basis1p)

cutoff(B1p::ACE.Product1pBasis) = cutoff(B1p.bases[2])

cutoff(Rn::ACE.Rn1pBasis) = cutoff(Rn.R)

# already defined in pair potentials 
# cutoff(Rn::ACE.OrthPolys.TransformedPolys) = Rn.ru


"""
`struct ACESitePotential` : wraps an ACE model into a JuLIP calculator
"""
struct ACESitePotential{ENV, TM} <: SitePotential
   models::Dict{AtomicNumber, TM}
end

struct ACESitePotentialBasis{ENV, TM} <: JuLIP.MLIPs.IPBasis
   models::Dict{AtomicNumber, TM}  # model = basis
   inds::Dict{AtomicNumber, UnitRange{Int}}
end

ACESiteCalc{ENV} = Union{ACESitePotential{ENV}, ACESitePotentialBasis{ENV}}

Base.length(ipbasis::ACESitePotentialBasis) = sum(length, values(ipbasis.models))

cutoff(V::ACESiteCalc) =  maximum(cutoff, values(V.models))

function ACESitePotential(models::Dict{Symbol, TM}, 
                          ENV = AtomicEnvironment) where {TM} 
   mods = Dict([ AtomicNumber(sym) => mo for (sym, mo) in models]...)
   return ACESitePotential(mods, ENV)
end

ACESitePotential(models::Dict{AtomicNumber, TM}, ENV = AtomicEnvironment) where {TM} = 
     ACESitePotential{ENV, TM}(models)

function _get_basisinds(V::ACESitePotential)
   inds = Dict{AtomicNumber, UnitRange{Int}}()
   i0 = 0
   for (z, mo) in V.models
      len = length(mo.basis)
      inds[z] = (i0+1):(i0+len)   # to generalize for general models
      i0 += len
   end
   return inds 
end

function basis(V::ACESitePotential{ENV}) where ENV 
   models = Dict( [sym => model.basis for (sym, model) in V.models]... )
   inds = _get_basisinds(V)
   return ACESitePotentialBasis{ENV, valtype(models)}(models, inds)
end

function environment(V::ACESiteCalc{<: AtomicEnvironment}, 
                     at::Atoms, nlist, i::Integer) where {TM}
   Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, i)
   z0 = at.Z[i]
   return environment(V, Rs, Zs, z0), Js
end

atomstate(z::AtomicNumber) = AtomState{Float64}( (mu = z,) )
atomstate(z::AtomicNumber, rr::SVector{3, T}) where {T} = 
      AtomState{T}( (mu = z, rr = rr) )

environment(V::ACESiteCalc{<: AtomicEnvironment}, 
            Rs::AbstractVector, Zs::AbstractVector, z0::AtomicNumber) = 
      AtomicEnvironment( atomstate(z0), atomstate.(Zs, Rs) )


JuLIP.alloc_temp(V::ACESiteCalc, N::Integer) = 
      ( JuLIP.Potentials.alloc_temp_site(N)..., 
      )

evaluate!(tmp, V::ACESitePotential, Rs, Zs, z0) = 
      evaluate!(tmp.tmpmodel[z0], V.models[z0], environment(V, Rs, Zs, z0)).val

function evaluate!(B, tmp, V::ACESitePotentialBasis, Rs, Zs, z0) 
   # fill!(B, 0)
   Bview = (@view B[V.inds[z0]])
   evaluate!(Bview, tmp.tmpmodel[z0], V.models[z0], environment(V, Rs, Zs, z0))
   return B 
end

# alloc_temp_d(V::ACESiteCalc, env::AbstractConfiguration) = 
#       ( JuLIP.Potentials.alloc_temp_site(length(env))...,
#         dV = zeros(JVec{Float64}, length(env)), 
#         tmpdmodel = Dict([ z => alloc_temp_d(mo, env) for (z, mo) in V.models ]...)
#       )

import ACEbase
function ACEbase.evaluate_d(V::ACESiteCalc, Rs::AbstractVector{JVec{T}}, Zs, z0) where {T} 
   env = environment(V, Rs, Zs, z0)
   tmpd = alloc_temp_d(V, env)
   ACE.grad_config!(tmpd.dV, tmpd.tmpdmodel[z0], V.models[z0], env)
end

function evaluate_d!(dV, _tmpd, V::ACESitePotential, Rs, Zs, z0) 
   env = environment(V, Rs, Zs, z0)
   tmpd = alloc_temp_d(V.models[z0], env)
   return ACE.grad_config!(dV, tmpd, V.models[z0], env)
end

function evaluate_d!(dB, _tmpd, V::ACESitePotentialBasis, Rs, Zs, z0)
   fill!(dB, zero(eltype(dB)))
   dBview = (@view dB[V.inds[z0], 1:length(Rs)])
   env = environment(V, Rs, Zs, z0)
   tmpd = alloc_temp_d(V.models[z0], env)
   evaluate_d!(dBview, tmpd, V.models[z0], env)
   return dBview
end
                
                

ACE.nparams(V::ACESitePotential) = sum(ACE.nparams, values(V.models))

function ACE.params(V::ACESitePotential) 
   inds = _get_basisinds(V)
   c = zeros(ACE.nparams(V))
   for (z, mo) in V.models 
      c[inds[z]] .= ACE.params(mo)
   end
   return c 
end 

function ACE.set_params!(V::ACESitePotential, c::AbstractVector)
   inds = _get_basisinds(V)
   for (z, mo) in V.models
      ACE.set_params!(V.models[z], c[inds[z]])
   end
   return V
end