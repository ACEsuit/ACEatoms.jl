@testset "ACEatoms.jl" begin

##

using ACE, ACEatoms, JuLIP, LinearAlgebra, ACEbase, Test
using ACEatoms: Species1PBasis, ZμRnYlm_1pbasis, AtomState,
                rand_environment, AtomicEnvironment
using ACEbase.Testing
using ACE: evaluate, evaluate_d, alloc_temp, alloc_B
using StaticArrays

inc_env(env, Us) = AtomicEnvironment(env.X0, env.Xs .+ Us)
notdot(a, b) = sum(notdot.(a, b))
notdot(a::Number, b::Number) = a * b

##

@info("Basic evaluation test for debugging")
Nat = 2
B1p = ACEatoms.ZμRnYlm_1pbasis(; species = [:C,:O])
env = rand_environment(B1p, Nat)
A1 = ACE.evaluate(B1p, env)
A, dA = ACE.evaluate_ed(B1p, env)
println(@test(A ≈ A1))
println(@test(size(dA) == (length(A), Nat)))
println(@test(eltype(dA) <: SVector))

##

@info("Range of finite difference tests")
for species in (:X, :Si, [:Ti, :Al], [:C, :H, :O])
   @info("    species $(species)")
   B1p = ACEatoms.ZμRnYlm_1pbasis(; species = species)
   # test deserialization
   # TODO Testing.test_fio

   for ntest = 1:20
      Nat = rand(5:15)
      env = rand_environment(B1p, Nat)
      Us = randn(SVector{3, Float64}, Nat)
      c = randn(length(B1p))
      F = t -> notdot(c, ACE.evaluate(B1p, inc_env(env, t * Us)))
      dF = t -> notdot(sum(Diagonal(c) * ACE.evaluate_ed(B1p, inc_env(env, t*Us))[2], dims = (1,))[:], Us)
      print_tf(@test(all(ACEbase.Testing.fdtest(F, dF, 0.0; verbose=false))))
   end
   println()
end


##
end
