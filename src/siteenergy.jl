

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

function basis(V::ACESitePotential{ENV}) where ENV 
   models = Dict( [sym => model.basis for (sym, model) in V.models]... )
   inds = Dict{AtomicNumber, UnitRange{Int}}()
   i0 = 0
   for (z, mo) in models
      inds[z] = (i0+1):(i0+length(mo))
      i0 += length(mo)
   end
   return ACESitePotentialBasis{ENV, valtype(models)}(models, inds)
end

function environment(V::ACESiteCalc{<: AtomicEnvironment}, 
                     at::Atoms, nlist, i::Integer) where {TM}
   Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, i)
   z0 = at.Z[i]
   return environment(V, Rs, Zs, z0), Js
end

environment(V::ACESiteCalc{<: AtomicEnvironment}, 
            Rs::AbstractVector, Zs::AbstractVector, z0::AtomicNumber) = 
      AtomicEnvironment( AtomState(z0), AtomState.(Zs, Rs) )


JuLIP.alloc_temp(V::ACESiteCalc, N::Integer) = 
      ( JuLIP.Potentials.alloc_temp_site(N)..., 
        tmpmodel = Dict([z => alloc_temp(mo) for (z, mo) in V.models]...)
      )

evaluate!(tmp, V::ACESitePotential, Rs, Zs, z0) = 
      evaluate!(tmp.tmpmodel[z0], V.models[z0], environment(V, Rs, Zs, z0)).val

function evaluate!(B, tmp, V::ACESitePotentialBasis, Rs, Zs, z0) 
   # fill!(B, 0)
   Bview = (@view B[V.inds[z0]])
   evaluate!(Bview, tmp.tmpmodel[z0], V.models[z0], environment(V, Rs, Zs, z0))
   return B 
end


alloc_temp_d(V::ACESiteCalc, N::Integer) = 
      ( JuLIP.Potentials.alloc_temp_site(N)...,
        dV = zeros(JVec{Float64}, N), 
        tmpdmodel = Dict([ z => alloc_temp_d(mo, N) for (z, mo) in V.models ]...)
      )

evaluate_d!(dV, tmpd, V::ACESitePotential, Rs, Zs, z0) = 
      ACE.grad_config!(dV, tmpd.tmpdmodel[z0], V.models[z0], environment(V, Rs, Zs, z0))

function evaluate_d!(dB, tmpd, V::ACESitePotentialBasis, Rs, Zs, z0)
   fill!(dB, zero(eltype(dB)))
   dBview = (@view dB[V.inds[z0], :])
   evaluate_d!(dBview, tmpd.tmpdmodel[z0], V.models[z0], environment(V, Rs, Zs, z0))
   return dB
end
                
                


