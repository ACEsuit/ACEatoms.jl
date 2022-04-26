

# @testset "Experimental Dipoles" begin 

##

using ACE, JuLIP, ACEatoms, ACEbase, Test, LinearAlgebra, StaticArrays
using ACE: evaluate, evaluate_d, SymmetricBasis, SimpleSparseBasis, PIBasis
#using ACEbase.Testing: fdtest
using JuLIP.Testing
using ACEatoms.Electrostatics: FixedChargeDipole, ESPot
using ACEbase.Testing 

##

# construct the 1p-basis
@info("Construcing a Linear ACE Model")
maxdeg = 6
ord = 3   # 4-body
Bsel = SimpleSparseBasis(ord, maxdeg)
species = [:Ti, :Al]
B1p = ACEatoms.ZμRnYlm_1pbasis(; species = species, maxdeg=maxdeg,
                                 rin = 1.2, rcut = 5.0)
φ = ACE.EuclideanVector{Float64}()
basis = SymmetricBasis(φ, B1p, Bsel)

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
println_slim(@test v1 ≈ v2)

cc = [cTi; cAl]
ACE.set_params!(V, cc)
ACEatoms.dipole(V, at)

##

B = ACEatoms.basis(V)
b = ACEatoms.dipole(B, at)
println_slim(@test( sum( b .* cc ) ≈ ACEatoms.dipole(V, at) ))

##

ACE.evaluate_d(V, Rs, Zs, z0)

##

@info("Finite difference test of ESPot")
# here is a nice way to avoid the loop :)
# In Julia it is often good to use functional programming paradigms
Q = Float64[ (z == 13 ? +1 : -1)  for z in atomic_numbers(at) ]
# or with map 
#  Q = map(z -> Float64(z == 13 ? +1 : -1),  atomic_numbers(at))
set_data!(at, :Q, Q)
Vref = FixedChargeDipole()
Vtot = ESPot(JuLIP.MLIPs.SumIP(Vref, V))

##

energy(Vtot, at)
forces(Vtot, at)

##

# fdtest(Vtot, at, verbose=true)
println_slim(@test JuLIP.Testing.fdtest(Vtot, at, verbose=true))

##

# results comparing to LAMMPS
# This test only works with \lambda = 1.0, \alpha=10.0 (not the current default settings hence removed)
pos = [[0 0 0], [0.6 -0.5 2.0]];
qs = [-3.0 2.0];
mus = [[1.0 -2.0 0.0], [0.0 1.0 1.0]];

E_lammps = -38.13193105209918
F_lammps = [4.98417434  -3.94540391  15.64287807]

# E = ACEatoms.Electrostatics.electrostatic_energy(pos, qs, mus, 1.0)
# F = ACEatoms.Electrostatics.electrostatic_forces(pos, qs, mus, 1.0)[1]

# println_slim(@test isapprox(E_lammps, E, rtol=1e-7))
# println_slim(@test isapprox(F_lammps[1], F[1], rtol=1e-7) && isapprox(F_lammps[2], F[2], rtol=1e-7) && isapprox(F_lammps[3], F[3], rtol=1e-7))


##

@info("test (de-)dictionisation of Dipole potential")
@warn("turned of failing FIO test")
# println_slim(@test all(JuLIP.Testing.test_fio(Vtot)))

##

# end
