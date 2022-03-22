



##

using ACE, JuLIP, ACEatoms, ACEbase, Test, LinearAlgebra
using ACE: evaluate, evaluate_d, SymmetricBasis, SimpleSparseBasis, PIBasis
using ACEbase.Testing: fdtest, println_slim 


##

# construct the 1p-basis
@info("Construcing a Linear ACE Model")
maxdeg = 8
ord = 3
Bsel = SimpleSparseBasis(ord, maxdeg)

species = [:Ti, :Al]
B1p = ACEatoms.ZμRnYlm_1pbasis(; species = species, maxdeg=maxdeg, Bsel = Bsel, 
                                 rin = 1.2, rcut = 5.0)

ACE.init1pspec!(B1p, Bsel)
basis = SymmetricBasis(ACE.Invariant(), B1p, Bsel)
cTi = randn(length(basis))
cAl = randn(length(basis))
models = Dict(:Ti => ACE.LinearACEModel(basis, cTi; evaluator = :standard), 
              :Al => ACE.LinearACEModel(basis, cAl; evaluator = :standard) )

@info("Convert LinearACEModel into an ACESitePotential")
V = ACEatoms.ACESitePotential(models)

## FIO 

@info("Check FIO")
using ACEbase.Testing: test_fio 
println_slim(@test(all(test_fio(V; warntype = false))))

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
println_slim(@test v1 ≈ v2)
energy(V, at)


##
@info("Test that grad_config(model) ≈ evaluate_d(site-potential)")
dEs = zeros(JVecF, length(env))
dv1 = ACE.grad_config!(dEs, V.models[z0], env)
dv2 = JuLIP.evaluate_d!(dEs, nothing, V, Rs, Zs, z0)
dv3 = evaluate_d(V, Rs, Zs, z0)

println_slim(@test(dv1 ≈ dv2))
println_slim(@test(dv1 ≈ dv3))

forces(V, at)
##

@info("Finite-difference test on total energy")
println_slim(@test JuLIP.Testing.fdtest(V, at))

##

@info("Check also the ACESitePotentialBasis interface")

ipbasis = ACEatoms.basis(V)
# convention is that z_Al < z_Ti hence this ordering 
cc = [cAl; cTi]

# check get_params, set_params

@info("  ... energy")
val1 = energy(V, at)
val2 = sum(cc .* energy(ipbasis, at))
println_slim(@test (val1 ≈ val2))

@info("  ... forces")
val1 = forces(V, at)
frcB = forces(ipbasis, at)
val2 = sum(cc .* frcB)
println_slim(@test val1 ≈ val2)

@info("  ... virial")
val1 = virial(V, at)
virB = virial(ipbasis, at)
val2 = sum(cc .* virB)
println_slim(@test val1 ≈ val2)
