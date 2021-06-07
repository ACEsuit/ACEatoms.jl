
using JuLIP: neighbourlist, cutoff
using ACEatoms.Electrostatics: Fixed_q_dipole, c_light, e
using JuLIP

_myreal(mu::ACE.EuclideanVector) = real.(mu.val)

_myreal(mu::AbstractVector{<: ACE.EuclideanVector}) = _myreal.(mu)

# experimental dipole function 
# this assumes that the underlying ACE model is covariant and 
# won't perform any checks
function dipole(V::ACESitePotential, at::Atoms)
   nlist = neighbourlist(at, cutoff(V))
   mu = sum( evaluate(V.models[at.Z[i]], 
                        environment(V, at, nlist, i)[1])
             for i = 1:length(at) )
   return _myreal(mu)
end


function dipole(B::ACESitePotentialBasis, at::Atoms)
   nlist = neighbourlist(at, cutoff(B))
   MU = zeros(SVector{3, Float64}, length(B))
   for i = 1:length(at)
      z0 = at.Z[i]
      env, _ = environment(B, at, nlist, i)
      mu = evaluate(B.models[z0], env) 
      MU[B.inds[z0]] .+= _myreal(mu)
   end
   return MU
end

"""
JuLIP dipole calculator returning the dipole moment of the atomic charges
"""
function dipole(Vref::Fixed_q_dipole, at::Atoms)
   mu = zeros(SVector{3, Float64})
   if has_data(at, :Q)
      Q = get_data(at, :Q)::Vector{Float64}
      mu += sum(Q .* positions(at) .* (1e-11/c_light/e), dims = 1)[1]
   end
   if has_data(at, :mu)
      at_mus = get_data(at, :mu)::Vector{SVector{3, Float64}}
      mu += sum(at_mus, dims = 1)[1]
   end
   return mu
end

function dipole(IP::JuLIP.MLIPs.SumIP{Any}, at::Atoms{Float64})
   mu = zeros(SVector{3, Float64})
   for pot in IP.components
      mu += dipole(pot, at)
   end
   return mu
end

function atomic_dipole!(V::ACESitePotential, at::Atoms)
   mu = zeros(SVector{3, Float64}, length(at) )
   nlist = neighbourlist(at, cutoff(V))
   for i = 1:length(at)
      mu[i] = _myreal(evaluate(V.models[at.Z[i]], 
                        environment(V, at, nlist, i)[1]))
   end
   set_data!(at, :mu, mu)
end
