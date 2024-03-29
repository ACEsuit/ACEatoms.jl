
##

using ACE, ACEatoms, JuLIP, LinearAlgebra, ACEbase, Test
using ACEatoms: Species1PBasis, ZμRnYlm_1pbasis, AtomState,
                rand_ACEConfig, PopZμRnYlm_1pbasis, rand_ACEConfig_pop
using ACEbase.Testing
using ACE: evaluate, evaluate_d
using StaticArrays



inc_env(env, Us) = ACEConfig(env.Xs .+ Us)
notdot(a, b) = sum(notdot.(a, b))
notdot(a::Number, b::Number) = a * b

##

@info("Basic evaluation test for debugging")
Nat = 2
B1p = ACEatoms.ZμRnYlm_1pbasis(; maxdeg = 10, species = [:C,:O])
env = rand_ACEConfig(B1p, Nat)
A1 = ACE.evaluate(B1p, env)
A, dA = ACE.evaluate_ed(B1p, env)
println_slim(@test(A ≈ A1))
println_slim(@test(size(dA) == (length(A), Nat)))

##

@info("Basic evaluation test for debugging of Pop basis set")
maxdeg = 6
Bsel = ACE.SimpleSparseBasis(1, maxdeg)
B1p_pop = ACEatoms.PopZμRnYlm_1pbasis(; maxdeg=maxdeg, species = [:C,:O])
B1p = ACE.Product1pBasis( tuple(B1p_pop.bases[2:end]...) )
ACE.init1pspec!(B1p, Bsel)
env = rand_ACEConfig_pop(B1p_pop, Nat)
A_pop = ACE.evaluate(B1p_pop, env)
A = ACE.evaluate(B1p, env)
println_slim(@test length(A) == length(A_pop))

# manual implementation: 
A1 = sum( evaluate(B1p, X) * X.pop for X in env.Xs )
println_slim(@test(A1 ≈ A_pop))

##
@info("Range of finite difference tests")
for species in (:X, :Si, [:Ti, :Al], [:C, :H, :O])
   @info("    species $(species)")
   local B1p = ACEatoms.ZμRnYlm_1pbasis(; species = species, maxdeg = 10, )
   # test deserialization
   # TODO Testing.test_fio
   _rrval(x) = x.rr

   for ntest = 1:20
      local Nat, env, F
      Nat = rand(5:15)
      env = rand_ACEConfig(B1p, Nat)
      Us = randn(SVector{3, Float64}, Nat)
      c = randn(length(B1p))
      F = t -> notdot(c, ACE.evaluate(B1p, inc_env(env, t * Us)))
      dF = t -> notdot( _rrval.(sum(Diagonal(c) * ACE.evaluate_d(B1p, inc_env(env, t*Us)), dims = (1,))[:]), Us)
      print_tf(@test(all(ACEbase.Testing.fdtest(F, dF, 0.0; verbose=false))))
   end
   println()
   B1p = ACEatoms.PopZμRnYlm_1pbasis(; species = species, maxdeg = 8)
   # test deserialization
   # TODO Testing.test_fio

   for ntest = 1:20
      local Nat, env, F
      Nat = rand(5:15)
      env = rand_ACEConfig_pop(B1p, Nat)
      Us = randn(SVector{3, Float64}, Nat)
      c = randn(length(B1p))
      F = t -> notdot(c, ACE.evaluate(B1p, inc_env(env, t * Us)))
      dF = t -> notdot( _rrval.(sum(Diagonal(c) * ACE.evaluate_d(B1p, inc_env(env, t*Us)), dims = (1,))[:]), Us)
      print_tf(@test(all(ACEbase.Testing.fdtest(F, dF, 0.0; verbose=false))))
   end
   println()
end

##
