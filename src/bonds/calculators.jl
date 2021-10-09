


function prepare_bondcalc(model, at)
   rcut_b = cutoff(model)
   rcut_e = cutoff_env(model)
   nlist_b = neighbourlist(at, rcut_b)
   nlist_e = neighbourlist(at, rcut_e)
   # consider filtering the nlist_b so that each bond 
   # occurs only ones, i.e. duplicates are deleted 
   return nlist_b, nlist_e
end

function get_bond_env(i, rri, j, rrj, nlist_e, at)
   Jsi, Rsi = neigs(nlist_e, i)
   Jsj, Rsj = neigs(nlist_e, j)
   Js = Int[] 
   Xs = JVecF[] 
   # go through all, fill them into a NearestNeighbour datastructure 
   rrm = 0.5 * (rri + rrj)
   for (Js0, Rs0, rr0) in ( (Jsi, Rsi, rri), (Jsj, Rsj, rrj) )
      for (j, rr) in zip(Js0, Rs0)
         # recenter the rr vector relative to bond mid-point
         rr = rr + rr0 - rrm
         # todo: if already in list skip 
         # todo: add it to the list 
         # convert it to an ACEState! add other stuff, such as species, etc. ???
         X = State(rr = rr, be = :e)
         if filter(bondmodel, X
            push!(Js, j)
            push!(Xs, X)
         end 
      end 
   end
   return Js, ACEConfig(Xs)
end



function energy(model::ACEBondPotential, at)
   nlist_b, nlist_e = prepare_bondcalc(model, at)

   for i = 1:length(at) 
      Js, Rs = neigs(nlist, i)
      for (j, rrij) in Rs 
         # todo: we want something like this but it is not obvious 
         # how to implement it 
         # if j > i; continue; end 

         # figure out what rri, rrj means in this context!?!?! not obvious 
         # again!!! 
         Js_env, cfg_env = get_bond_env(i, rri, j, rrj, nlist_e, at)

         # add the center-bond to the config 
         push!(cfg.Xs, State(rr = rrij, be = b))
         
         # evaluate the ace bond model and add to total energy 
         E += evaluate(model, cfg) 
      end 
   end 

   return E 
end

