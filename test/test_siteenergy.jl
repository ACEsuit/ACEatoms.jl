

using ACE, JuLIP, ACEatoms, ACEbase, Test, LinearAlgebra
using ACE: evaluate, evaluate_d, SymmetricBasis, NaiveTotalDegree, PIBasis
using ACEbase.Testing: fdtest


##

# construct the 1p-basis
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

V = ACEatoms.ACESitePotential(model)

##

zTi = AtomicNumber(:Ti)
zAl = AtomicNumber(:Al)
at = bulk(:Ti, cubic=true) * 3
at.Z[2:3:end] .= zAl
at = rattle!(at, 0.1)

##

nlist = neighbourlist(at, cutoff(V))
env, _ = ACEatoms.environment(V, at, nlist, 2)
Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, 2)
z0 = at.Z[2]
env1 = ACEatoms.environment(V, Rs, Zs, z0)
ACEbase._allfieldsequal(env1, env)
evaluate(model, env)
evaluate(V, Rs, Zs, z0)

##

JuLIP.usethreads!(false)
JuLIP.energy(V, at)
