

using ACE, JuLIP, ACEatoms, ACEbase, Test, LinearAlgebra
using ACE: evaluate, evaluate_d, SymmetricBasis, NaiveTotalDegree, PIBasis
using ACEbase.Testing: fdtest

# I2LModel = ACEatoms.Models.I2LModel

##

# construct the 1p-basis
@info("Construcing L1")
D = NaiveTotalDegree()
maxdeg = 5
ord = 3
species = [:Ti, :Al]
B1p = ACEatoms.ZμRnYlm_1pbasis(; species = species, maxdeg=maxdeg, D = D, 
                                 rin = 1.2, rcut = 5.0)
ACE.init1pspec!(B1p, maxdeg = maxdeg, Deg = ACE.NaiveTotalDegree())
φ = ACE.Invariant()
pibasis = PIBasis(B1p, ord, maxdeg; property = φ)
basis = SymmetricBasis(pibasis, φ)
cTi = randn(length(basis))
cAl = randn(length(basis))
L1 = Dict(:Ti => ACE.LinearACEModel(basis, cTi; evaluator = :standard), 
          :Al => ACE.LinearACEModel(basis, cAl; evaluator = :standard) )

##
@info("Construcing L2")
Pk_in = ACE.transformed_jacobi(maxdeg, ACE.Transforms.IdTransform(), 
                               1.0, 0.0; pin = 0, pcut = 0)
Pk = ACE.Scal1pBasis(Pk_in, :u, :k)
Bnlmk = B1p * Pk
ACE.init1pspec!(Bnlmk, maxdeg = maxdeg, Deg = ACE.NaiveTotalDegree())
length(Bnlmk)
φ = ACE.Invariant()
pibasis2 = PIBasis(Bnlmk, ord, maxdeg; property = φ)
basis2 = SymmetricBasis(pibasis2, φ)
cTi = randn(length(basis2))
cAl = randn(length(basis2))
L2 = Dict(:Ti => ACE.LinearACEModel(basis2, cTi; evaluator = :standard), 
          :Al => ACE.LinearACEModel(basis2, cAl; evaluator = :standard) )

##

@info("Construct I2LModel")
model = ACEatoms.Models.I2LModel(L1, L2)

@info("Basic parameter tests")
println(@test( ACE.nparams(model) == 2 * length(basis) + 2 * length(basis2) ))
c2 = randn(ACE.nparams(model))
ACE.set_params!(model, c2)
println(@test( ACE.params(model) == c2 ))

##

@info("construct random atomic structure")

zTi = AtomicNumber(:Ti)
zAl = AtomicNumber(:Al)
at = bulk(:Ti, cubic=true) * 3
at.Z[2:3:end] .= zAl
at = rattle!(at, 0.1)

##

Us = rand(length(at))
nlist = neighbourlist(at, 5.0)
env1, _ = ACEatoms.environment(model.L1, at, nlist, 1)
env2 = ACEatoms.Models.AtEnvU1(at, nlist, 1, Us)
ACE.evaluate(model.L1.models[at.Z[1]], env1)
ACE.evaluate(model.L2.models[at.Z[1]], env2)

##

energy(model, at)

## forces 


## derivatives w.r.t. parameters 

