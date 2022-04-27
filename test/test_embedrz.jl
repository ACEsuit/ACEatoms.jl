
using ACE, ACEatoms, JuLIP, LinearAlgebra, ACEbase, Test
using ACEatoms: AtomState
using ACEbase.Testing
using ACE: evaluate, evaluate_d, State
using ACE: @λ
using StaticArrays
using ACEatoms: embed_nl_spec1, embedrz

##

# create a radial basis 
maxdeg = 10
r0 = 1.0
rcut = 3.0
trans = @λ r -> 2/(1+r)
Qn = ACE.transformed_jacobi(maxdeg, trans, rcut)
# this is a bad temporary hack, I need to figure out a way to get this 
# more easily ... it's all happening, just needs a bit of time.
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
evaluate_d(Rnl, X)

## basic tests 

@info("check correct application of the embedding, correctness of `evaluate`")
for ntest = 1:30 
   X = rand_X()
   print_tf(@test(embeddings[X.mu] * evaluate(Qn, norm(X.rr)) ≈ evaluate(Rnl, X)))
end

##

using ACEbase.Testing: fdtest 
@info("Finite difference test")
for ntest = 1:30 
   u = randn(length(Rnl))
   dX = ACE.DState( rr = ACE.rand_sphere() )
   F = t -> dot(u, evaluate(Rnl, X + t * dX))
   dF = t -> (g = evaluate_d(Rnl, X + t * dX);  dot(dX.rr, sum(u[n] * g[n].rr for n = 1:length(g)))) 
   print_tf(@test( fdtest(F, dF, 0.0; verbose=false) ))
end

##

