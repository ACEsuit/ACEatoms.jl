

##

using ACEatoms, Test 
using ACEbase.Testing
ES = ACEatoms.Electrostatics

##

@info("Total dipole moment test")
X = [0.0 0.0 0.0; 1.0 0.0 0.0]
Q = [-1.0; 1.0]
MU = [0.0 0.0 0.0; 0.0 1.0 -2.0]
tot_dipole = ES.get_dipole(X, Q, MU)
println_slim(@test(isapprox(tot_dipole, [0.208194 1.0 -2.0], atol=1e-6) ))

##

@info("Soft charge-charge interaction energy test")
X = [0.0 0.0; 1.0 0.0]
q1 = +1.0
q2 = -1.0
E = ES.soft_coulomb(X[1, :], q1, X[2, :], q2, 1.0)
println_slim(@test( isapprox(E, -14.39964547, atol=1e-6) )) 

##

@info("Soft charge-dipole interaction energy test")
X = [0.0 0.0; 1.0 0.0]
q1 = 1.0
dip = [1.0 0.0]
E = ES.soft_q_μ(X[1, :], q1, X[2, :], dip, 1.0)
println_slim(@test( isapprox(E, -2.99792458, atol=1e-6) ))

##

