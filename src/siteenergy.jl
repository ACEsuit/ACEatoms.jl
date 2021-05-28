

import JuLIP.Potentials: SitePotential

# TODO: nasty hack - must fix this!!!

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
   model::TM
end

struct ACESitePotentialBasis{ENV, TM} <: JuLIP.MLIPs.IPBasis
   model::TM  # model = basis
end

ACESiteCalc{ENV} = Union{ACESitePotential{ENV}, ACESitePotentialBasis{ENV}}


Base.length(ipbasis::ACESitePotentialBasis) = length(ipbasis.model)

cutoff(V::ACESitePotential) =  cutoff(V.model.basis)
cutoff(V::ACESitePotentialBasis) =  cutoff(V.model)


ACESitePotential(model, ENV = AtomicEnvironment) = 
     ACESitePotential{ENV, typeof(model)}(model)

basis(V::ACESitePotential, ENV = AtomicEnvironment) = 
      ACESitePotentialBasis{ENV, typeof(V.model.basis)}(V.model.basis)

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
        tmpmodel = alloc_temp(V.model)
      )

evaluate!(tmp, V::ACESitePotential, Rs, Zs, z0) = 
      evaluate!(tmp.tmpmodel, V.model, environment(V, Rs, Zs, z0)).val

evaluate!(B, tmp, V::ACESitePotentialBasis, Rs, Zs, z0) = 
      evaluate!(B, tmp.tmpmodel, V.model, environment(V, Rs, Zs, z0))
      

alloc_temp_d(V::ACESiteCalc, N::Integer) = 
      ( JuLIP.Potentials.alloc_temp_site(N)...,
        dV = zeros(JVec{Float64}, N), 
        tmpdmodel = alloc_temp_d(V.model, N)
      )

evaluate_d!(dV, tmpd, V::ACESitePotential, Rs, Zs, z0) = 
      ACE.grad_config!(dV, tmpd.tmpdmodel, V.model, environment(V, Rs, Zs, z0))

evaluate_d!(dB, tmpd, V::ACESitePotentialBasis, Rs, Zs, z0) = 
      evaluate_d!(dB, tmpd.tmpdmodel, V.model, environment(V, Rs, Zs, z0))
                
                

# ----------------- similar interface for the basis 


