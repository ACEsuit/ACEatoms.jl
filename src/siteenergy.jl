

import JuLIP.Potentials: SitePotential

function cutoff(B1p::ACE.Product1pBasis)
   Rn = B1p.bases[2]
   return Rn.R.ru
end

"""
`struct ACESitePotential` : wraps an ACE model into a JuLIP calculator
"""
struct ACESitePotential{ENV, TM} <: SitePotential
   model::TM
end

cutoff(V::ACESitePotential) =  cutoff(V.model.basis.pibasis.basis1p)


ACESitePotential(model, ENV = AtomicEnvironment) = 
     ACESitePotential{ENV, typeof(model)}(model)


function environment(V::ACESitePotential{<: AtomicEnvironment}, 
                     at::Atoms, nlist, i::Integer) where {TM}
   Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, i)
   z0 = at.Z[i]
   return AtomicEnvironment( AtomState(z0), AtomState.(Zs, Rs) ), Js
end

environment(V::ACESitePotential{<: AtomicEnvironment}, 
            Rs::AbstractVector, Zs::AbstractVector, z0::AtomicNumber) = 
      AtomicEnvironment( AtomState(z0), 
                         AtomState.(Zs, Rs) )

JuLIP.alloc_temp(V::ACESitePotential, N::Integer) = 
      ( JuLIP.Potentials.alloc_temp_site(N)..., 
        tmpmodel = alloc_temp(V.model)
      )


evaluate!(tmp, V::ACESitePotential, Rs, Zs, z0) = 
      evaluate!(tmp.tmpmodel, V.model, environment(V, Rs, Zs, z0)).val

      
alloc_temp_d(V::ACESitePotential, N::Integer) = 
      ( JuLIP.Potentials.alloc_temp_site(N)...,
        dV = zeros(JVec{Float64}, N), 
        tmpdmodel = alloc_temp_d(V.model, N)
      )

evaluate_d!(dV, tmpd, V::ACESitePotential, Rs, Zs, z0) = 
      ACE.grad_config!(dV, tmpd.tmpdmodel, V.model, 
                      environment(V, Rs, Zs, z0))



