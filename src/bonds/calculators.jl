
import ACE: cutoff_env 

function prepare_bondcalc(model, at)
   rcut_b = cutoff(model)
   rcut_e = cutoff_env(model)    #  CAREFUL HOW WE DEFINE THIS!!!!!
   nlist_b = neighbourlist(at, rcut_b)
   nlist_e = neighbourlist(at, rcut_e)
   # consider filtering the nlist_b so that each bond 
   # occurs only ones, i.e. duplicates are deleted 
   return nlist_b, nlist_e
end

function get_bond_env(i0, rr0, i1, rr1, nlist_e, bondmodel)
   Js0, Rs0 = neigs(nlist_e, i0)  # superset of actual environment
   Js = Int[] 
   Xs = JVecF[] 
   # go through all, fill them into a NearestNeighbour datastructure 
   # FIX THIS : 0.5 -> lambda in bondmodel
   rrm = 0.5 * (rr0 + rr1)
   for (j, rr_ij) in zip(Js0, Rs0)
      # recenter the rr vector relative to bond mid-point
      rr_j = rr0 + rr_ij   # absolute position of that env atom 
      rr_mj = rr_j - rrm   # position relative  to bond center point 
      # convert it to an ACEState add other stuff, such as species, etc. ???
      X = State(rr = rr_mj, be = :e, rr0 = rr1 - rr0)
      if filter(bondmodel, X)
         push!(Js, j)
         push!(Xs, X)
      end 
   end

   # add the center-bond to the config 
   push!(Js, 0)
   push!(cfg.Xs, State(rr = rr1 - rr0, be = :b, rr0 = rr1 - rr0))

   return Js, ACEConfig(Xs)
end



function energy(model::ACEBondPotential, at)
   nlist_b, nlist_e = prepare_bondcalc(model, at)

   for i = 1:length(at) 
      Js, Rs = neigs(nlist_b, i)
      for (j, rr_ij) in Rs 
         # todo: we want something like this but it is not obvious 
         # how to implement it : if j > i; continue; end 

         # figure out what rri, rrj means in this context!?!?! not obvious 
         # again!!! 
         rr_i = at.X[i] 
         rr_j = rr_i + rr_ij  # relative position; j could be outside the cell
         Js_env, cfg_env = get_bond_env(i, rr_i, j, rr_j, nlist_e, bondmodel)
         
         # evaluate the ace bond model and add to total energy 
         E += evaluate(model, cfg) 
      end 
   end 

   return E 
end

