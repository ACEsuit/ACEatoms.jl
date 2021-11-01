"""
`module Electrostatics`

This module implements various functions related to (long range) electrostatic interactions
"""

module Electrostatics

using LinearAlgebra
using Zygote
using JuLIP
using StaticArrays
using ACEatoms, ACE
using SpecialFunctions: erf
import JuLIP: energy, forces
import ACE: write_dict, read_dict

export get_dipole, electrostatic_energy, electrostatic_forces, soft_coulomb, 
        soft_q_μ, soft_μ_μ, soft_q_μ2, dipole, FixedChargeDipole, energy, forces,
        write_dict, read_dict

const μ_0 = 4.0e-7 * π
const c = 299792458.0
const ϵ_0 = 1 / μ_0 / c^2
const e = 1.602176634e-19
const c_light = 299792458

"""
  JuLIP calculator using a dipole ACE model to calculate the soft core Coulomb Q-Q, Q-μ, μ-μ 
  energies and forces
"""

struct ESPot{TEV} <: AbstractCalculator 
  dipoleevaluator::TEV
end

write_dict(V::ESPot) = Dict("__id__" => "ACEatoms_ESPot",
                           "dipoleevaluator" => write_dict(V.dipoleevaluator))

function read_dict(::Val{:ACEatoms_ESPot}, D::Dict) 
  dipoleevaluator = read_dict(D["dipoleevaluator"])
  return ESPot(dipoleevaluator)
end

==(V1::ESPot, V2::ESPot) = 
      all(V1.dipoleevaluator.components .== V2.dipoleevaluator.components)

function energy(V::ESPot, at::Atoms) 
  # get dipoles with dipole evaluator, and than get energy
  pos = positions(at)
  mus = ACEatoms.atomic_dipole(V.dipoleevaluator.components[2], at)
  if has_data(at, :Q)
    Qs = get_data(at, :Q)::Vector{Float64} 
  end
  return electrostatic_energy(pos, Qs, mus, 0.0)  # Set λ=0.0 as a default, may need to change later
end

function forces(V::ESPot, at::Atoms) 
  λ = 0.0  # a sensible looking default for the soft core parameter
  pos = positions(at)
  # CO: I don't like `V.dipoleevaluator.components[2]` at all, it feels very 
  #     fragile. -> to discuss!!        
  nlist = neighbourlist(at, cutoff(V.dipoleevaluator.components[2]))
  # TODO getting the atomic dipole at both the enenergy and force calls is wasteful if one calls them together
  #      CO: I don't bother with that in the Julia code. The evaluation is usually
  #          a fraction (~20-30%) of the gradients. Yes, we can optimise this 
  #          but only worth for production...
  mus = ACEatoms.atomic_dipole(V.dipoleevaluator.components[2], at)
  
  # get the static charges on each atom 
  if has_data(at, :Q)
    Qs = get_data(at, :Q)::Vector{Float64} 
  else 
    Qs = zeros(length(at))  
  end

  # allocate vectors for storing forces, evaluating the potentials, etc
  Fs = zeros(JVec{Float64}, length(at))   # forces 
  fs = zeros(JVec{Float64}, length(at))   # fs[i] <- ∂E / ∂μ[i]
  maxN = JuLIP.maxneigs(nlist)
  # tmpd = ACE.alloc_temp_d(V.dipoleevaluator.components[2], maxN)
  # dVc = zeros(SMatrix{3,3,ComplexF64}, maxN)
  dV = zeros(SMatrix{3,3,Float64}, maxN)

  for i = 1:length(at)
    for j = (i+1):length(at) 
      Rij = pos[i] - pos[j]
      fqμ1 = force_q_μ(Rij, Qs[i], mus[j], λ)  # -> (DRij, Dmu)
      Fs[i] -= fqμ1[1]
      Fs[j] += fqμ1[1]
      fs[j] += fqμ1[2]
      fqμ2 = force_q_μ(-1.0 .* Rij, Qs[j], mus[i], λ)
      Fs[i] += fqμ2[1]
      Fs[j] -= fqμ2[1]
      fs[i] += fqμ2[2]
      fμμ = force_μ_μ(Rij, mus[i], mus[j], λ)
      Fs[i] -= fμμ[1]
      Fs[j] += fμμ[1]  
      fs[i] += fμμ[2]
      fs[j] += fμμ[3]
      fqq = force_q_q(Rij, Qs[i], Qs[j], λ)[1]
      Fs[i] -= fqq
      Fs[j] += fqq
    end 
  end
  
  for i = 1:length(at)
    Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, i)
    z0 = at.Z[i]
    # ACE.grad_config(dV, V.dipoleevaluator.components[2], Rs, Zs, z0)
    # map!(dv -> real.(dv), dV, dVc)
    dV = ACE.evaluate_d(V.dipoleevaluator.components[2], Rs, Zs, z0)

    for (j, dVij) in zip(Js, dV)
      Fs[j] -= transpose(dVij) * fs[i]
      Fs[i] += transpose(dVij) * fs[i]
    end
  end 
  return Fs
end


"""
Get the overall dipole moment of a set of point charges and point dipoles
Arguments:
  pos::Array natoms x 3 Array
  charges::Array natoms x 1 Array
  dipoles::Array natoms x 3 Array
Returns:
  total_dipole::Array 1 x 3 Array
"""
function get_dipole(pos::AbstractArray, charges::AbstractVector, dipoles::AbstractArray, pbc::Bool=false)
  @assert pbc == false "Periodic boundary condition not yet supported"
  return sum(charges .* pos .* (1e-11/c_light/e) .+ dipoles, dims=1)
end

"""
Fixed charge (or dipole) reference potential
"""
struct FixedChargeDipole end

write_dict(::FixedChargeDipole) = Dict("__id__" => "ACEatoms_FixedChargeDipole")

read_dict(::Val{:ACEatoms_FixedChargeDipole}, D::Dict) = FixedChargeDipole()

"""
Total electrostatic enenrgy of a set of point charges and point dipoles
calculated using the soft core potentials. 
λ = 1.0 returns the non-soft core version.
"""
function electrostatic_energy(pos::AbstractArray, charges::AbstractArray, dipoles::AbstractArray, λ::Real, pbc::Bool=false)
  @assert (!pbc) "Periodic boundary condition not yet supported"
  qq = 0.0
  qμ = 0.0
  μμ = 0.0
  for (i, R) in enumerate(pos)
    for j = (i+1):length(pos)
      Rij = R - pos[j]
      qq += soft_coulomb(Rij, charges[i], charges[j], λ)
      qμ += soft_q_μ(Rij, charges[i], dipoles[j], λ)
      qμ += soft_q_μ(-1.0 .* Rij, charges[j], dipoles[i], λ)
      μμ += soft_μ_μ(Rij, dipoles[i], dipoles[j], λ)
    end
  end
  return qq + qμ + μμ
end

"""
Electrostatic forces of a set of FIXED point charges and FIXED point dipoles
calculated using the soft core potentials. 
λ = 1.0 returns the non-soft core version.
"""
function electrostatic_forces(pos::AbstractArray, charges::AbstractArray, dipoles::AbstractArray, λ::Real, pbc::Bool=false)
  if pbc; error("Periodic boundary condition not yet supported"); end 
  Fs = zeros(JVec{Float64}, length(charges))
  for (i, R) in enumerate(pos)
    for j = (i+1):length(pos)
      Rij = R - pos[j]
      fqq = force_q_q(Rij,charges[i], charges[j], λ)[1]
      Fs[i] -= JVec(fqq)
      Fs[j] += JVec(fqq)
      fqμ = force_q_μ(Rij, charges[i], dipoles[j], λ)[1]
      Fs[i] -= JVec(fqμ)
      Fs[j] += JVec(fqμ)
      fqμ = force_q_μ(-1.0 .* Rij, charges[j], dipoles[i], λ)[1]
      Fs[i] += JVec(fqμ)
      Fs[j] -= JVec(fqμ)
      fμμ = force_μ_μ(Rij, dipoles[i], dipoles[j], λ)[1]
      Fs[i] -= JVec(fμμ)
      Fs[j] += JVec(fμμ)    
    end
  end
  return Fs
end

"""
Soft core Coulomb interaction between two point charges
If λ=1.0 the normal Coulomb is returned
Adapted from LAMMPS `pair_style coul/long/soft`
"""
function soft_coulomb(Rij::AbstractArray, q1::Real, q2::Real, λ::Real=0.0, α::Real=10.0)
  return q1 * q2 / (4*π*ϵ_0 * (α * (1 - λ)^2 + dot(Rij, Rij))^0.5) * e * 1e10
end

function force_q_q(Rij::AbstractArray, q1::Real, q2::Real, λ::Real=0.0, α::Real=10.0)
  return Zygote.gradient(r -> soft_coulomb(r, q1, q2, λ, α), Rij)
end

"""
Coulomb interaction between two Gaussian charge clouds
"""

function gaus_coulomb(Rij::AbstractArray, q1::Real, q2::Real, σ1::Real, σ2::Real)
  rij = dot(Rij, Rij)^0.5
  α = 1 / (2 * (σ1^2 + σ2^2))^0.5
  return q1 * q2 / (4*π*ϵ_0 * rij) * erf(α*rij) * e * 1e10
end

function force_q_q_gauss(Rij::AbstractArray, q1::Real, q2::Real, σ1::Real, σ2::Real)
  return Zygote.gradient(r -> gaus_coulomb(r, q1, q2, σ1, σ2), Rij)
end

"""
Soft core Coulomb interaction between point charge and point dipole
If λ=1.0 the normal Coulomb is returned
Analogously defined to LAMMPS `pair_style coul/long/soft`
Rji = Ri - Rj vector
"""
function soft_q_μ(Rji::AbstractArray, q1::Real, μ::AbstractArray, λ::Real=0.0, α::Real=10.0)
  return q1 * dot(μ, Rji) / (4*π*ϵ_0 * (α * (1 - λ)^2 + dot(Rji, Rji))^1.5) * (1e-1/c_light)
end

function force_q_μ(Rji::AbstractArray, q1::Real, μ::AbstractArray, λ::Real=0.0, α::Real=10.0)
  return Zygote.gradient((r, mu) -> soft_q_μ(r, q1, mu, λ, α), Rji, μ)
end

"""
Soft core Coulomb interaction between two point dipoles
If λ=1.0 the normal Coulomb is returned
Analogously defined to LAMMPS `pair_style coul/long/soft`
Rji = Ri - Rj vector
"""
function soft_μ_μ(Rji::AbstractArray, μ1::AbstractArray, μ2::AbstractArray, λ::Real=0.0, α::Real=10.0)
  return 1.0 / (4*π*ϵ_0 * (α * (1.0 - λ)^2 + dot(Rji, Rji))^1.5) * (dot(μ1, μ2) - 3 * dot(μ1, Rji) * dot(μ2, Rji) / (α * (1 - λ)^2 + dot(Rji, Rji))) *1e-12 * (1/c_light)^2 / e
end

function force_μ_μ(Rji::AbstractArray, μ1::AbstractArray, μ2::AbstractArray, λ::Real=0.0, α::Real=10.0)
  return Zygote.gradient((r, mu1, mu2) -> soft_μ_μ(r, mu1, mu2, λ, α), Rji, μ1, μ2)
end

end # end of module Electrostatics