
import ACE
export PolyPairBasis

import ACEbase: valtype, acquire_B!, acquire_dB!, 
                  release_B!, release_dB!


# TODO: allow PairBasis with different Pr for each z, z' combination

function pairbasis(species, maxdeg::Integer, rcut, trans; pcut = 2)
   Pr = transformed_jacobi(maxdeg, trans, rcut; pcut = pcut)
   # Pr is now a chain ...
   len = length(Pr.F[end])
   zlist = ZList(species; static=true)
   bidx0 = get_bidx0(Pr, zlist)
   return PolyPairBasis(Pr, zlist, bidx0, len, rcut)
end

struct PolyPairBasis{TJ, NZ} <: IPBasis
   J::TJ
   zlist::SZList{NZ}
   bidx0::SMatrix{NZ,NZ,Int}
   len::Int
   rcut::Float64 
end

valtype(pB::PolyPairBasis, args...) = valtype(pB.J)

Base.length(pB::PolyPairBasis) = pB.len * (numz(pB) * (numz(pB) + 1)) ÷ 2
Base.length(pB::PolyPairBasis, z0::AtomicNumber) = length(pB, z2i(pB, z0))
Base.length(pB::PolyPairBasis, iz0::Integer) = pB.len

zlist(pB::PolyPairBasis) = pB.zlist

function scaling(pB::PolyPairBasis, p)
   ww = zeros(Float64, length(pB))
   for iz0 = 1:numz(pB), iz = 1:numz(pB)
      idx0 = _Bidx0(pB, iz0, iz)
      for n = 1:length(pB)
         # TODO: very crude, can we do better?
         #       -> need a proper H2-orthogonbality?
         ww[idx0+n] = n^p
      end
   end
   return ww
end

PolyPairBasis(J, species, len, rcut) = 
   PolyPairBasis( J, ZList(species; static=true), len, rcut)

PolyPairBasis(J, zlist::SZList, rcut) =
   PolyPairBasis(J, zlist, get_bidx0(J, zlist), len, rcut)

function get_bidx0(J, zlist::SZList{NZ}) where {NZ}
   NJ = length(J)
   bidx0 = fill(zero(Int), (NZ, NZ))
   i0 = 0
   for i = 1:NZ, j = i:NZ
      bidx0[i,j] = i0
      bidx0[j,i] = i0
      i0 += NJ
   end
   return SMatrix{NZ, NZ, Int}(bidx0...)
end


==(B1::PolyPairBasis, B2::PolyPairBasis) =
      (B1.J == B2.J) && (B1.zlist == B2.zlist) && (B1.rcut == B2.rcut)

JuLIP.cutoff(pB::PolyPairBasis) = pB.rcut 

write_dict(pB::PolyPairBasis) = Dict(
      "__id__" => "ACE_PolyPairBasis",
          "Pr" => write_dict(pB.J),
       "zlist" => write_dict(pB.zlist) )

read_dict(::Val{:ACE_PolyPairBasis}, D::Dict) =
      PolyPairBasis( read_dict(D["Pr"]), read_dict(D["zlist"]) )

# alloc_temp(pB::PolyPairBasis, args...) = (
#               J = alloc_B(pB.J),
#           tmp_J = alloc_temp(pB.J)  )

# alloc_temp_d(pB::PolyPairBasis, args...) =  (
#              J = alloc_B( pB.J),
#          tmp_J = alloc_temp(pB.J),
#             dJ = alloc_dB(pB.J),
#         tmpd_J = alloc_temp_d(pB.J)  )

"""
compute the zeroth index of the basis corresponding to the potential between
two species zi, zj; as precomputed in `PolyPairBasis.bidx0`
"""
_Bidx0(pB, zi, zj) = pB.bidx0[ z2i(pB, zi), z2i(pB, zj) ]
_Bidx0(pB, i::Integer, j::Integer) = pB.bidx0[ i, j ]

function energy(pB::PolyPairBasis, at::Atoms{T}) where {T}
   E = zeros(T, length(pB))
   for (i, j, R) in pairs(at, cutoff(pB))
      r = norm(R)
      J = evaluate(pB.J, r)
      idx0 = _Bidx0(pB, at.Z[i], at.Z[j])
      for n = 1:length(pB)
         E[idx0 + n] += 0.5 * J[n]
      end
      ACE.release!(J)
   end
   return E
end

function forces(pB::PolyPairBasis, at::Atoms{T}) where {T}
   F = zeros(JVec{T}, length(at), length(pB))
   for (i, j, R) in pairs(at, cutoff(pB))
      r = norm(R)
      dJ = evaluate_d(pB.J, r)
      idx0 = _Bidx0(pB, at.Z[i], at.Z[j])
      for n = 1:length(pB)
         F[i, idx0 + n] += 0.5 * dJ[n] * (R/r)
         F[j, idx0 + n] -= 0.5 * dJ[n] * (R/r)
      end
      ACE.release!(dJ)
   end
   return [ F[:, iB] for iB = 1:length(pB) ]
end

function virial(pB::PolyPairBasis, at::Atoms{T}) where {T}
   V = zeros(JMat{T}, length(pB))
   for (i, j, R) in pairs(at, cutoff(pB))
      r = norm(R)
      dJ = evaluate_d(pB.J, r)
      idx0 = _Bidx0(pB, at.Z[i], at.Z[j])
      for n = 1:length(pB)
         V[idx0 + n] -= 0.5 * (dJ[n]/r) * R * R'
      end
      ACE.release!(dJ)
   end
   return V
end
