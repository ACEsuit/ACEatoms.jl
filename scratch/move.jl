

# TODO: move to ACEatoms
# # -> rand_config?
# function rand_nhd(Nat, J::ScalarBasis, species = :X)
#    zlist = ZList(species)
#    Rs = [ rand_vec(J) for _ = 1:Nat ]
#    Zs = [ rand(zlist.list) for _ = 1:Nat ]
#    z0 = rand(zlist.list)
#    return Rs, Zs, z0
# end

# rand_config(species; kwargs...) =
#       rand_config(ZList(species); kwargs...)

# rand_config(V::AbstractCalculator; kwargs...) =
#       rand_config(zlist(V); kwargs...)

# function rand_config(zlist::Union{ZList, SZList};
#                      absrattle = 0.0, relrattle = 0.2, repeat = 3,
#                      kwargs...)
#    # start with the longest rnn
#    rnns = rnn.(zlist.list)
#    sym = chemical_symbol( zlist.list[findmax(rnns)[2]] )
#    at = bulk(sym; kwargs...) * repeat
#    for n = 1:length(at)
#       at.Z[n] = rand(zlist.list)
#    end
#    rattle!(at, maximum(rnns) * relrattle + absrattle)
#    return at
# end



# @info("Check several properties of PIBasis")
# for species in (:X, :Si, [:C, :O, :H]), N = 1:5
#    local AA, AAdag, dAA, dAAdag, dagbasis, basis, Rs, Zs, z0
#    maxdeg = 7
#    Nat = 15
#    P1 = ACE.RnYlm1pBasis(Pr; species = species)
#    dagbasis = ACE.PIBasis(P1, N, D, maxdeg)
#    basis = standardevaluator(dagbasis)
#    @info("species = $species; N = $N; length = $(length(basis))")
#    @info("test (de-)serialisation")
#    # only require dagbasis to deserialize correctly
#    println(@test all(JuLIP.Testing.test_fio(dagbasis)))
#    @info("Check Permutation invariance")
#    for ntest = 1:20
#       Rs, Zs, z0 = ACE.rand_nhd(Nat, Pr, species)
#       p = randperm(length(Rs))
#       print_tf(@test(evaluate(basis, Rs, Zs, z0) ≈
#                      evaluate(basis, Rs[p], Zs[p], z0)))
#    end
#    println()
#    @info("Check gradients")
#    for ntest = 1:20
#       Rs, Zs, z0 = ACE.rand_nhd(Nat, Pr, species)
#       AA = evaluate(basis, Rs, Zs, z0)
#       dAA = evaluate_d(basis, Rs, Zs, z0)
#       Us = [ rand(eltype(Rs)) .- 0.5 for _=1:length(Rs) ]
#       dAA_dUs = transpose.(dAA) * Us
#       errs = []
#       for p = 2:12
#          h = 0.1^p
#          AA_h = evaluate(basis, Rs + h * Us, Zs, z0)
#          dAA_h = (AA_h - AA) / h
#          # @show norm(dAA_h - dAA_dUs, Inf)
#          push!(errs, norm(dAA_h - dAA_dUs, Inf))
#       end
#       success = (/(extrema(errs)...) < 1e-3) || (minimum(errs) < 1e-10)
#       print_tf(@test success)
#    end
#
#    println()
#    @info("Check Standard=DAG Evaluator")
#    for ntest = 1:20
#       Rs, Zs, z0 = ACE.rand_nhd(Nat, Pr, species)
#       AA = evaluate(basis, Rs, Zs, z0)
#       AAdag = evaluate(dagbasis, Rs, Zs, z0)
#       print_tf(@test AA ≈ AAdag)
#
#       dAA = evaluate_d(basis, Rs, Zs, z0)
#       dAAdag = evaluate_d(dagbasis, Rs, Zs, z0)
#       print_tf(@test dAA ≈ dAAdag)
#    end
#    println()
# end
# println()

#---





# #---
# @info("Basis construction and evaluation checks")
# @info("check single species")
# Nat = 15
# Rs, Zs, z0 = rand_nhd(Nat, Pr, :X)
# B = evaluate(rpibasis, Rs, Zs, z0)
# println(@test(length(rpibasis) == length(B)))
# dB = evaluate_d(rpibasis, Rs, Zs, z0)
# println(@test(size(dB) == (length(rpibasis), length(Rs))))
# B_, dB_ = evaluate_ed(rpibasis, Rs, Zs, z0)
# println(@test (B_ ≈ B) && (dB_ ≈ dB))
#
# #---
# @info("check multi-species")
# maxdeg = 5
# Pr = transformed_jacobi(maxdeg, trans, rcut; pcut = 2)
# species = [:C, :O, :H]
# P1 = ACE.RnYlm1pBasis(Pr; species = species, D = D)
# basis = ACE.RPIBasis(P1, N, D, maxdeg)
# Rs, Zs, z0 = ACE.rand_nhd(Nat, Pr, species)
# B = evaluate(basis, Rs, Zs, z0)
# println(@test(length(basis) == length(B)))
# dB = evaluate_d(basis, Rs, Zs, z0)
# println(@test(size(dB) == (length(basis), length(Rs))))
# B_, dB_ = evaluate_ed(basis, Rs, Zs, z0)
# println(@test (B_ ≈ B) && (dB_ ≈ dB))
#
# #---
#
# degrees = [ 12, 10, 8, 8, 8, 8 ]
#
# @info("Check a few basis properties ")
# # for species in (:X, :Si) # , [:C, :O, :H])
# for species in (:X, :Si, [:C, :O, :H]), N = 1:length(degrees)
#    local Rs, Zs, z0, B, dB, basis, D, P1, Nat
#    Nat = 15
#    D = SparsePSHDegree()
#    P1 = ACE.RnYlm1pBasis(Pr; species = species)
#    basis = ACE.RPIBasis(P1, N, D, degrees[N])
#    @info("species = $species; N = $N; deg = $(degrees[N]); len = $(length(basis))")
#    @info("   check (de-)serialization")
#    println(@test(all(JuLIP.Testing.test_fio(basis))))
#    @info("   isometry and permutation invariance")
#    for ntest = 1:30
#       Rs, Zs, z0 = ACE.rand_nhd(Nat, Pr, species)
#       Rsp, Zsp = ACE.rand_sym(Rs, Zs)
#       print_tf(@test(evaluate(basis, Rs, Zs, z0) ≈
#                      evaluate(basis, Rsp, Zsp, z0)))
#    end
#    println()
#    @info("   check derivatives")
#    for ntest = 1:30
#       Rs, Zs, z0 = ACE.rand_nhd(Nat, Pr, species)
#       B = evaluate(basis, Rs, Zs, z0)
#       dB = evaluate_d(basis, Rs, Zs, z0)
#       Us = [ rand(eltype(Rs)) .- 0.5 for _=1:length(Rs) ]
#       dB_dUs = transpose.(dB) * Us
#       errs = []
#       for p = 2:12
#          h = 0.1^p
#          B_h = evaluate(basis, Rs + h * Us, Zs, z0)
#          dB_h = (B_h - B) / h
#          # @show norm(dAA_h - dAA_dUs, Inf)
#          push!(errs, norm(dB_h - dB_dUs, Inf))
#       end
#       success = (/(extrema(errs)...) < 1e-3) || (minimum(errs) < 1e-10)
#       print_tf(@test success)
#    end
#    println()
#    @info("   check combine")
#    coeffs = randcoeffs(basis)
#    V = combine(basis, coeffs)
#    Vst = standardevaluator(V)
#    for ntest = 1:30
#       Rs, Zs, z0 = ACE.rand_nhd(Nat, Pr, species)
#       v = evaluate(V, Rs, Zs, z0)
#       vst = evaluate(Vst, Rs, Zs, z0)
#       cdotB = dot(coeffs, evaluate(basis, Rs, Zs, z0))
#       print_tf(@test v ≈ cdotB ≈ vst)
#    end
#    println()
#    @info("   check graph evaluator")
#    basisst = standardevaluator(basis)
#    for ntest = 1:30
#       env = ACE.rand_nhd(Nat, Pr, species)
#       print_tf(@test evaluate(basisst, env...) ≈ evaluate(basis, env...))
#       print_tf(@test evaluate_d(basisst, env...) ≈ evaluate_d(basis, env...))
#    end
#    println()
# end
#
