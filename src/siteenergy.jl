

import JuLIP.Potentials: SitePotential
import ACE: ACEConfig 
import ACEbase: OneParticleBasis

# TODO: nasty hack - must fix this!!!

cutoff(model::ACE.LinearACEModel) = cutoff(model.basis)

cutoff(basis::ACE.SymmetricBasis) = cutoff(basis.pibasis)

cutoff(basis::ACE.PIBasis) = cutoff(basis.basis1p)

cutoff(B1p::ACE.Product1pBasis) = minimum(cutoff.(B1p.bases))

cutoff(B1p::OneParticleBasis) = Inf

cutoff(Rn::ACE.B1pComponent) = 
      haskey(Rn.meta, "rcut") ? Rn.meta["rcut"]::Float64 : Inf

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
                          ENV = ACEConfig) where {TM} 
   mods = Dict([ AtomicNumber(sym) => mo for (sym, mo) in models]...)
   return ACESitePotential(mods, ENV)
end

ACESitePotential(models::Dict{AtomicNumber, TM}, ENV = ACEConfig) where {TM} = 
     ACESitePotential{ENV, TM}(models)



# ----------------------------- FIO 
import ACEbase: write_dict, read_dict
import Base: == 

write_dict(V::ACESitePotential) = 
      Dict( "__id__" => "ACEatoms_ACESitePotential", 
            "models" => Dict([ string(Int(z)) => write_dict(m)  for (z, m) in V.models ]...))

function read_dict(::Val{:ACEatoms_ACESitePotential}, D::Dict) 
   models = Dict([ AtomicNumber(parse(Int, z)) => read_dict(m) 
                   for (z, m) in D["models"] ]...)
   return ACESitePotential(models)
end

==(V1::ACESitePotential, V2::ACESitePotential) = 
      (V1.models == V2.models)

# ----------------------------- 


function _get_basisinds(V::ACESitePotential)
   inds = Dict{AtomicNumber, UnitRange{Int}}()
   zz = sort(collect(keys(V.models)))
   i0 = 0
   for z in zz
      mo = V.models[z]
      len = length(mo.basis)
      inds[z] = (i0+1):(i0+len)   # to generalize for general models
      i0 += len
   end
   return inds 
end

_get_basisinds(V::ACEatoms.ACESitePotentialBasis) = V.inds


function basis(V::ACESitePotential{ENV}) where ENV 
   models = Dict( [sym => model.basis for (sym, model) in V.models]... )
   inds = _get_basisinds(V)
   return ACESitePotentialBasis{ENV, Base.valtype(models)}(models, inds)
end

function environment(V::ACESiteCalc{<: ACEConfig}, 
                     at::Atoms, nlist, i::Integer) where {TM}
   Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, i)
   z0 = at.Z[i]
   return environment(V, Rs, Zs, z0), Js
end

atomstate(z::AtomicNumber) = AtomState{Float64}( (mu0 = z,) )
atomstate(z::AtomicNumber, z0::AtomicNumber, rr::SVector{3, T}) where {T} = 
      AtomState{T}( (mu = z, mu0 = z0, rr = rr) )

environment(V::ACESiteCalc{<: ACEConfig}, 
            Rs::AbstractVector, Zs::AbstractVector, z0::AtomicNumber) = 
      ACEConfig( atomstate.(Zs, z0, Rs) )


evaluate!(tmp, V::ACESitePotential, Rs, Zs, z0) = 
      evaluate(V.models[z0], environment(V, Rs, Zs, z0)).val

# JuLIP.alloc_temp(V::ACESitePotentialBasis, N::Integer) = 
#       (JuLIP.Potentials.alloc_temp_site(N)..., )

function evaluate!(B, tmp, V::ACESitePotentialBasis, Rs, Zs, z0) 
   # fill!(B, 0)
   Bview = (@view B[V.inds[z0]])
   evaluate!(Bview, V.models[z0], environment(V, Rs, Zs, z0))
   return B 
end

import ACEbase
function ACEbase.evaluate_d(V::ACESiteCalc, Rs::AbstractVector{JVec{T}}, Zs, z0) where {T} 
   env = environment(V, Rs, Zs, z0)
   dV = ACE.grad_config(V.models[z0], env)
   return [ dv.rr for dv in dV ] 
end

function evaluate_d!(dV, _tmpd, V::ACESitePotential, Rs, Zs, z0) 
   env = environment(V, Rs, Zs, z0)
   g = ACE.grad_config(V.models[z0], env)
   for i = 1:length(g)
      dV[i] = g[i].rr
   end
   ACE.release!(g)
   return dV
end

function evaluate_d!(dB, _tmpd, V::ACESitePotentialBasis, Rs, Zs, z0)
   fill!(dB, zero(eltype(dB)))
   dBview = (@view dB[V.inds[z0], 1:length(Rs)])
   env = environment(V, Rs, Zs, z0)
   g = evaluate_d(V.models[z0], env) 
   for i = 1:length(g)
      dBview[i] = g[i].rr
   end
   return dBview
end

                

ACE.nparams(V::ACESitePotential) = sum(ACE.nparams, values(V.models))
ACE.nparams(V::ACEatoms.ACESitePotentialBasis) = length(V)

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

function ACE.scaling(V::ACESitePotentialBasis, p)
   inds = _get_basisinds(V)
   scal = Vector{Float64}(undef, ACE.nparams(V))
   for (z, model) in V.models
      scal[inds[z]] .= ACE.scaling(model, p)
   end
   return scal
end