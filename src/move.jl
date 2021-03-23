

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
