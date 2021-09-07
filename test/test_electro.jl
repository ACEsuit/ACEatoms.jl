@testset "Electrostatics" begin

##

using ACEatoms, Test 
ES = ACEatoms.Electrostatics

##

@info("Total dipole moment test")
positions = [0.0 0.0 0.0; 1.0 0.0 0.0]
charges = [-1.0; 1.0]
dipoles = [0.0 0.0 0.0; 0.0 1.0 -2.0]
tot_dipole = ES.get_dipole(positions, charges, dipoles)
println(@test(isapprox(tot_dipole, [0.208194 1.0 -2.0], atol=1e-6) ))

##

@info("Soft charge-charge interaction energy test")
positions = [0.0 0.0; 1.0 0.0]
q1 = +1.0
q2 = -1.0
E = ES.soft_coulomb(positions[1, :], q1, positions[2, :], q2, 1.0)
println(@test( isapprox(E, -14.39964547, atol=1e-6) )) 

##

@info("Soft charge-dipole interaction energy test")
positions = [0.0 0.0; 1.0 0.0]
q1 = 1.0
dip = [1.0 0.0]
E = ES.soft_q_μ(positions[1, :], q1, positions[2, :], dip, 1.0)
println(@test( isapprox(E, -2.99792458, atol=1e-6) ))

##


end