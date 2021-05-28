

@testset "ACESitePotential" begin 

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
c = randn(length(basis))
model = ACE.LinearACEModel(basis, c; evaluator = :standard)

@info("Convert LinearACEModel into an ACESitePotential")
V = ACEatoms.ACESitePotential(model)

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
env1 = ACEatoms.environment(V, Rs, Zs, z0)
ACEbase._allfieldsequal(env1, env)
v1 = evaluate(model, env).val
v2 = evaluate(V, Rs, Zs, z0)
println(@test v1 ≈ v2)
energy(V, at)

##
@info("Test that grad_config(model) ≈ evaluate_d(site-potential)")
dEs = zeros(JVecF, length(env))
dv1 = ACE.grad_config!(dEs, ACE.alloc_temp_d(V.model, env), V.model, env)
dv2 = JuLIP.evaluate_d!(dEs, ACE.alloc_temp_d(V, length(env)), V, Rs, Zs, z0)
dv3 = evaluate_d(V, Rs, Zs, z0)

println(@test(dv1 ≈ dv2))
println(@test(dv1 ≈ dv3))

forces(V, at)
##

@info("Finite-difference test on total energy")
println(@test JuLIP.Testing.fdtest(V, at))


##

end