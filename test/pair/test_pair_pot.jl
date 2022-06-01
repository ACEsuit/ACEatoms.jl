



@info("--------------- PolyPairPot Implementation ---------------")

##

using ACE, ACEatoms, ACEbase, Printf, Test, LinearAlgebra, JuLIP, JuLIP.Testing
using JuLIP: evaluate, evaluate_d
using JuLIP.Potentials: i2z, numz
using JuLIP.MLIPs: combine
using ACEatoms.PairPotentials: pairbasis
using ACEbase.Testing: fdtest 

randr() = 1.0 + rand()
randcoeffs(B) = rand(length(B)) .* (1:length(B)).^(-2)
JuLIP._usethreads[] = false 

##

maxdeg = 8
r0 = 1.0
rcut = 3.0

trans = polytransform(1, r0)
pB = pairbasis(:W, maxdeg, rcut, trans) 

coeffs = randcoeffs(pB)
V = combine(pB, coeffs)

## 

@info("FD test for pair potential in r space")
for ntest = 1:30 
   local coeffs, V 
   coeffs = randcoeffs(pB)
   V = combine(pB, coeffs)
   print_tf(@test fdtest(r -> ACE.evaluate(V, r), 
                         r -> ACE.evaluate_d(V, r), 
                         r0 + rand() * (rcut - r0); 
                         verbose=false))
end
println() 

##

at = bulk(:W, cubic=true) * 3
rattle!(at, 0.03)
r0 = rnn(:W)
X = copy(positions(at))
energy(V, at)

@info("Testing correctness of `PolyPairPot` against `PolyPairBasis`")
@info("    test `combine`")
coeffs = randcoeffs(pB)
V = combine(pB, coeffs)
_frcerr(F1, F2) = maximum(norm.(F1 - F2))
println(@test energy(V, at) ≈ sum(V.coeffs .*  energy(pB, at)))
println(@test _frcerr(forces(V, at), sum(coeffs .* forces(pB, at))) < 1e-12)

##

# frcs = [] 
# frcsb = [] 
# for n = 1:10 
#    push!(frcsb, sum(coeffs .* forces(pB, at)))
#    push!(frcs, forces(V, at))
# end

# ##
# all( frcs[n] == frcs[1] for n = 2:length(frcs) )
# frc_errs = [ _frcerr(frcs[n], frcs[1]) for n = 2:length(frcs) ]
# @show frc_errs

# # [ _frcerr(frcsb[n], frcsb[1]) for n = 2:length(frcs) ]

# ##
# vv = [] 
# for n = 1:10 
#    push!(vv, virial(V, at))
# end

# ##
# all( frcs[n] == frcs[1] for n = 2:length(frcs) )
# errs = [ norm(vv[n] - vv[1]) for n = 2:length(frcs) ]
# @show errs 


##

@info("   test (de-)dictionisation")
println(@test all(JuLIP.Testing.test_fio(V; warntype=false)))

@info("      check that PolyPairBasis ≈ PolyPairPot")
for ntest = 1:10
   local coeffs, V 
   rattle!(at, 0.01)
   coeffs = randcoeffs(pB)
   V = combine(pB, coeffs)

   E_V = energy(V, at)
   E_b = dot(energy(pB, at), coeffs)
   print_tf(@test E_V ≈ E_b)

   F_V = forces(V, at)
   F_b = sum(coeffs .* forces(pB, at))
   print_tf(@test F_V ≈ F_b)

   V_V = virial(V, at)
   V_b = sum(coeffs .* virial(pB, at))
   print_tf(@test V_V ≈ V_b)
end
println()

##

@info("      Standard JuLIP Force Consistency Test")
variablecell!(at)
rattle!(at, 0.03)
println(@test JuLIP.Testing.fdtest(V, at))

##
