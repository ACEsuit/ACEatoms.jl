

@testset "ACESitePotential" begin 

##

using ACE, JuLIP, ACEatoms, ACEbase, Test, LinearAlgebra
using ACE: evaluate, evaluate_d, SymmetricBasis, NaiveTotalDegree, PIBasis
using ACEbase.Testing: fdtest


##

# construct the 1p-basis
@info("Construcing a Linear ACE Model")
D = NaiveTotalDegree()
maxdeg = 8
ord = 3
species = [:Ti, :Al]
B1p = ACEatoms.ZμRnYlm_1pbasis(; species = species, maxdeg=maxdeg, D = D, 
                                 rin = 1.2, rcut = 5.0)
ACE.init1pspec!(B1p)
φ = ACE.Invariant()
pibasis = PIBasis(B1p, ord, maxdeg; property = φ)
basis = SymmetricBasis(pibasis, φ)
cTi = randn(length(basis))
cAl = randn(length(basis))
models = Dict(:Ti => ACE.LinearACEModel(basis, cTi; evaluator = :standard), 
              :Al => ACE.LinearACEModel(basis, cAl; evaluator = :standard) )

@info("Convert LinearACEModel into an ACESitePotential")
V = ACEatoms.ACESitePotential(models)

##

@info("Create random TiAl configuration")
zTi = AtomicNumber(:Ti)
zAl = AtomicNumber(:Al)
at = bulk(:Ti, cubic=true) * 3
at.Z[2:3:end] .= zAl
at = rattle!(at, 0.1)

##

@info("Test that evaluate(model) ≈ evaluate(site-potential)")
nlist = neighbourlist(at, cutoff(V))
env, _ = ACEatoms.environment(V, at, nlist, 2)
Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, 2)
z0 = at.Z[2]
sym = chemical_symbol(z0)
env1 = ACEatoms.environment(V, Rs, Zs, z0)
ACEbase._allfieldsequal(env1, env)
v1 = evaluate(models[sym], env).val
v2 = evaluate(V, Rs, Zs, z0)
println(@test v1 ≈ v2)
energy(V, at)

##
@info("Test that grad_config(model) ≈ evaluate_d(site-potential)")
dEs = zeros(JVecF, length(env))
dv1 = ACE.grad_config!(dEs, ACE.alloc_temp_d(V.models[z0], env), V.models[z0], env)
dv2 = JuLIP.evaluate_d!(dEs, ACE.alloc_temp_d(V, length(env)), V, Rs, Zs, z0)
dv3 = evaluate_d(V, Rs, Zs, z0)

println(@test(dv1 ≈ dv2))
println(@test(dv1 ≈ dv3))

forces(V, at)
##

@info("Finite-difference test on total energy")
println(@test JuLIP.Testing.fdtest(V, at))

##

@info("Check also the ACESitePotentialBasis interface")

ipbasis = ACEatoms.basis(V)
cc = [cTi; cAl]

# check get_params, set_params

@info("  ... energy")
val1 = energy(V, at)
val2 = sum(cc .* energy(ipbasis, at))
println(@test (val1 ≈ val2))

@info("  ... forces")
val1 = forces(V, at)
frcB = forces(ipbasis, at)
val2 = sum(cc .* frcB)
println(@test val1 ≈ val2)

@info("  ... virial")
val1 = virial(V, at)
virB = virial(ipbasis, at)
val2 = sum(cc .* virB)
println(@test val1 ≈ val2)

##

end