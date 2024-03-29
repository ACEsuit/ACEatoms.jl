


@info("-------- Test PolyPairBasis Implementation ---------")

##

using ACE, ACEatoms
using Printf, Test, LinearAlgebra, JuLIP, JuLIP.Testing
using JuLIP: evaluate, evaluate_d
using JuLIP.Potentials: i2z, numz

using ACEatoms.PairPotentials: pairbasis


randr() = 1.0 + rand()
randcoeffs(B) = rand(length(B)) .* (1:length(B)).^(-2)

##

maxdeg = 8
r0 = 1.0
rcut = 3.0

trans = polytransform(1, r0)

Pr = transformed_jacobi(maxdeg, trans, rcut; pcut = 2)
pB = pairbasis(:W, maxdeg, rcut, trans)

@info("Scaling Test")
println(@test ACE.scaling(pB, 1) == 1:length(pB))
println(@test ACE.scaling(pB, 2) == (1:length(pB)).^2)

##

at = bulk(:W, cubic=true) * 3
rattle!(at, 0.03)
r0 = rnn(:W)
X = copy(positions(at))

E = energy(pB, at)
println(@test length(E) == length(pB))
F = forces(pB, at)
println(@test (length(F) == length(pB)) && all(length.(F) .==  length(at)))


##

@info("test (de-)dictionisation of PairBasis")
println(@test all(JuLIP.Testing.test_fio(pB; warntype=false)))

##

@info("Finite-difference test on PolyPairBasis forces")
for ntest = 1:30
   local E = energy(pB, at)
   DE = - forces(pB, at)
   U = [ (rand(JVecF) .- 0.5) for _=1:length(at) ]
   DExU = dot.(DE, Ref(U))
   errs = Float64[]
   for p = 2:10
      h = 0.1^p
      Eh = energy(pB, set_positions!(at, X+h*U))
      DEhxU = (Eh-E) / h
      push!(errs, norm(DExU - DEhxU, Inf))
   end
   success = (/(extrema(errs)...) < 1e-3) || (minimum(errs) < 1e-10)
   print_tf(@test success)
end
println()
##

