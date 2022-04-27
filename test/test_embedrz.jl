
using ACE, ACEatoms, JuLIP, LinearAlgebra, ACEbase, Test
using ACEatoms: AtomState
using ACEbase.Testing
using ACE: evaluate, evaluate_d, State
using StaticArrays

using ACEatoms: embed_nl_spec1, embedrz

##

# create a radial basis 
maxdeg = 10
r0 = 1.0
rcut = 3.0
trans = PolyTransform(1, r0)
Qn = ACE.transformed_jacobi(maxdeg, trans, rcut)
# this is a bad temporary hack, I need to figure out a way to get this 
# more easily ...
len_Qn = length(Qn.F[end])

# now embed it 
species = [:Fe, :C]
# this is good for testing but we may want to extract the list of all 
# possible (n, l) tuples from the ACE basis.
spec = ACEatoms.embed_nl_spec1(5, 1.0)
len_Rnl = length(spec)
embeddings = Dict(:Fe => randn(len_Rnl, len_Qn), 
                   :C => randn(len_Rnl, len_Qn))
Rnl = embedrz(Qn, spec, embeddings)


# We can create a state 
rand_X() = State( rr = ACE.rand_sphere() * (r0 + rand() * (rcut-r0)), 
                  mu = rand(species) )

X = rand_X()
evaluate(Rnl, X)