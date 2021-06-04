


import ACEatoms 
import ACEatoms: AtomState, ACESitePotential

const AtomStateU1{T} = ACE.State{(:mu, :rr, :u), 
                                 Tuple{AtomicNumber, SVector{3, T}, T}}

const AtEnvDef = AtomicEnvironment{AtomState{Float64}}
const AtEnvU1 = AtomicEnvironment{AtomStateU1{Float64}}

struct I2LModel{TL1, TL2} <: AbstractACEModel
   L1::TL1
   L2::TL2
end

# ------------ Model construction codes

I2LModel(L1::Dict, L2::Dict) = 
      I2LModel( ACESitePotential(L1), ACESitePotential(L2) )


# ------------ Parameter wrangling 

nparams(m::I2LModel) = nparams(m.L1) + nparams(m.L2)

params(m::I2LModel) = [ params(m.L1); params(m.L2) ]

function set_params!(m::I2LModel, c)
   n1 = nparams(m.L1)
   set_params!(m.L1, c[1:n1])
   set_params!(m.L2, c[(n1+1):end])
end


# ----------- evaluation codes 

function AtEnvU1(at, nlist, i, U)
   Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, i)
   z0 = at.Z[i]
   X0 = AtomStateU1{Float64}( (mu = z0, u = U[i]) )
   Xs = [ AtomStateU1{Float64}( (mu = Zs[t], rr = Rs[t], u = U[Js[t]]) ) 
          for t = 1:length(Js) ]
   return AtomicEnvironment(X0, Xs)
end

alloc_temp(m::I2LModel, args...) = 
      ( tmpL1 = alloc_temp(m.L1, args...), 
        tmpL2 = alloc_temp(m.L2, args...)
      )

function JuLIP.energy(m::I2LModel, at::Atoms)
   # evaluate the first layer 
   nlist = JuLIP.neighbourlist(at, cutoff(m.L1))
   U = zeros(Float64, length(at))
   for i = 1:length(at) 
      env, Js = environment(m.L1, at, nlist, i)
      U[i] = evaluate(m.L1.models[at.Z[i]], env)
   end

   # filter the U values to go between 0, 1 
   map!(u -> 1 / (1 + exp(u)), U, U)

   # evaluate the second layer 
   E = zero(fltype(m.L2))
   nlist = JuLIP.neighbourlist(at, cutoff(m.L2 ))  # could technically be different rcut 
   for i = 1:length(at) 
      envU = AtEnvU1(at, nlist, i, U)
      E += evaluate(m.L2.models[at.Z[i]], envU)
   end

   return E.val
end
