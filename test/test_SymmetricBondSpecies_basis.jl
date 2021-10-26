using ACE, ACEatoms
using ProgressMeter
r0cut = 2.0
rcut = 1.0
zcut = 2.0
env = ACE.EllipsoidBondEnvelope(r0cut, rcut, zcut;floppy=false, λ= .5)

maxorder = 3
Bsel = ACE.PNormSparseBasis(maxorder; p = 2, default_maxdeg = 3) 

@info("Test invariance")

bond_basis = ACE.Utils.SymmetricBond_basis(ACE.EuclideanVector(), env, Bsel; )
@show length(bond_basis)
offsite_basis = ACEatoms.SymmetricBondSpecies_basis(ACE.EuclideanVector(), env, Bsel; species = [:Al,:Ti] )
@show length(offsite_basis)
#@show length(basis_cov)
