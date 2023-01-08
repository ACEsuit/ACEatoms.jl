
using ACE, ACEatoms, JuLIP, LinearAlgebra, ACEbase, StaticArrays
using ACEatoms: AtomState, embed_nl_spec1, embedrz
using ACEbase.Testing
using ACE: evaluate, evaluate_d, State, @λ

##

# create a radial basis from which to construct the embedding 
maxdeg = 10
r0 = 1.0
rcut = 3.0
trans = @λ r -> 2/(1+r)  # this is a new way to specify arbtirary transforms
Qn = ACE.transformed_jacobi(maxdeg, trans, rcut)

# We need to manually get the length of Qn; I've made some changes in ACE.jl
# and somehow this is something that got lost on the way :) - I'll fix this 
# eventually and then also fix this tutorial 
len_Qn = length(Qn.F[end])


# now embed it to obtain the actual Rn basis. 
species = [:Fe, :C]

# first we need to somehow generate the specification of the embedding. If we 
# only want an n channel (Rn) rather than the more general Rnl then we can simply 
# do 
maxn = 15
spec = [ (n = n,) for n = 1:maxn ]

# next we generate the actual embeddings 
len_Rn = length(spec)
embeddings = Dict(:Fe => randn(len_Rn, len_Qn), 
                   :C => randn(len_Rn, len_Qn))
Rn = embedrz(Qn, spec, embeddings)

# We can create a state and evaluate it 
rand_X() = State( rr = ACE.rand_sphere() * (r0 + rand() * (rcut-r0)), 
                  mu = rand(species) )

X = rand_X()
evaluate(Rn, X)
evaluate_d(Rn, X)

# To see what this really achieves - i.e. this really explains what is implemented: 
embeddings[X.mu] * evaluate(Qn, norm(X.rr)) ≈ evaluate(Rn, X)

## Now we can generate an ACE basis 

maxL = 6
Ylm = ACE.Ylm1pBasis(maxL)
B1p = ACE.Product1pBasis( (Rn, Ylm) )
Bsel = ACE.SimpleSparseBasis(3, maxn)

ACE.init1pspec!(B1p, Bsel)

function filter_embed(bb::AbstractVector{<: NamedTuple})
   if length(bb) == 0; return true; end 
   return all(b.n == bb[1].n for b in bb)
end

basis = ACE.SymmetricBasis(ACE.Invariant(), B1p, Bsel; filterfun = filter_embed)

# for some reason this takes a really long time to generate, not sure what is 
# going on here ...

cfg = ACE.ACEConfig([ rand_X() for _=1:10] )
evaluate(basis, cfg)
evaluate_d(basis, cfg)

## Another thought: suppose we want an Rnl embedding rather than an Rn embedding
# this would mean that each l channel gets different embeddings. this would look 
# as follows: 

maxn = 20
spec_nl = embed_nl_spec1(maxn, maxn/maxL)
len_Rnl = length(spec_nl)
embed_nl = Dict(:Fe => randn(len_Rnl, len_Qn), 
                 :C => randn(len_Rnl, len_Qn))
Rnl = embedrz(Qn, spec_nl, embed_nl)

B1p_nl = ACE.Product1pBasis( (Rnl, Ylm) )
Bsel_nl = ACE.SimpleSparseBasis(4, 20)

basis_nl = ACE.SymmetricBasis(ACE.Invariant(), B1p_nl, Bsel_nl; filterfun = filter_embed)

length(basis_nl)
evaluate(basis_nl, cfg)
evaluate_d(basis_nl, cfg)

# should now be able to easily build ACE models from this ... 
