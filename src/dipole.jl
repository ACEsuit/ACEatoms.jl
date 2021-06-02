
using JuLIP: neighbourlist, cutoff

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
