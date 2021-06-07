

@testset "Experimental Dipoles" begin 

##

using ACE, JuLIP, ACEatoms, ACEbase, Test, LinearAlgebra
using ACE: evaluate, evaluate_d, SymmetricBasis, NaiveTotalDegree, PIBasis
using ACEbase.Testing: fdtest


##

# construct the 1p-basis
@info("Construcing a Linear ACE Model")
D = NaiveTotalDegree()
maxdeg = 6
ord = 3   # 4-body
species = [:Ti, :Al]
B1p = ACEatoms.ZμRnYlm_1pbasis(; species = species, maxdeg=maxdeg, D = D, 
                                 rin = 1.2, rcut = 5.0)
ACE.init1pspec!(B1p, maxdeg = maxdeg, Deg = ACE.NaiveTotalDegree())
φ = ACE.EuclideanVector{ComplexF64}()
pibasis = PIBasis(B1p, ord, maxdeg; property = φ, isreal = false)
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

cc = [cTi; cAl]
ACE.set_params!(V, cc)
ACEatoms.dipole(V, at)

##

B = ACEatoms.basis(V)
b = ACEatoms.dipole(B, at)
println(@test( sum( b .* cc ) ≈ ACEatoms.dipole(V, at) ))

end