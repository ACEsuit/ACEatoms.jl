
##

using ACE, ACEatoms, JuLIP, LinearAlgebra, ACEbase, Test
using ACEatoms: Species1PBasis, ZμRnYlm_1pbasis, AtomState,
                rand_ACEConfig, PopZμRnYlm_1pbasis, rand_ACEConfig_pop
using ACEbase.Testing
using ACE: evaluate, evaluate_d, State 
using StaticArrays

using ACEatoms: EmbedCat1pBasis


##

categories = [:O, :H, :C, :N] 
embeddings = Dict( :O => randn(3), :H => randn(3), :C => randn(3), :N => randn(3) )
Sk = EmbedCat1pBasis(categories, embeddings; varsym = :μ, idxsym = :k, label = "Sk")

##

for ntest = 1:20 
   local X 
   μ = rand(categories)
   X = State( μ = μ, rr = 3 * randn(3) )
   _S = evaluate(Sk, X)
   print_tf(@test( _S == embeddings[μ] ))
end

##

# now construct a product basis 
RnYlm = ACE.Utils.RnYlm_1pbasis()
B1p = RnYlm * Sk
println(B1p["Rn"])
println(B1p["Ylm"])
# println_slim(@test B1p["Rn"] isa ACE.Rn1pBasis )
# println_slim(@test B1p["Ylm"] isa ACE.Ylm1pBasis )
println_slim(@test B1p["Sk"] isa EmbedCat1pBasis )

Bsel = ACE.SimpleSparseBasis(3, 5)
ACE.init1pspec!(B1p, Bsel)
ACE.init1pspec!(RnYlm, Bsel)

println_slim(@test length(B1p) == 3 * length(RnYlm) )

## 

μ = rand(categories)
X = State( μ = μ, rr = 3 * randn(3) )

evaluate(B1p, X) 
evaluate(RnYlm, X)
evaluate(Sk, X)
