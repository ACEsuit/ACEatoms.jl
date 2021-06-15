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
import JuLIP: energy, forces

export get_dipole, electrostatic_energy, electrostatic_forces, soft_coulomb, 
        soft_q_μ, soft_μ_μ, soft_q_μ2, dipole, FixedChargeDipole, energy, forces

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
  nlist = neighbourlist(at, cutoff(V.dipoleevaluator.components[2]))
  # TODO getting the atomic dipole at both the enenergy and force calls is wasteful if oe calls them together
  mus = ACEatoms.atomic_dipole(V.dipoleevaluator.components[2], at)
  Qs = zeros(length(at))
  if has_data(at, :Q)
    Qs = get_data(at, :Q)::Vector{Float64} 
  end  
  Fs = zeros(JVec{Float64}, length(at))
  fs = zeros(JVec{Float64}, length(at))
  for i = 1:length(at)
    Js, Rs, Zs = JuLIP.Potentials.neigsz(nlist, at, i)
    z0 = at.Z[i]
    tmpd = ACE.alloc_temp_d(V.dipoleevaluator.components[2], length(Rs))
    dV = zeros(SMatrix{3,3,ComplexF64}, length(Rs))
    ACE.evaluate_d!(dV, tmpd, V.dipoleevaluator.components[2], Rs, Zs, z0)
    for j = (i+1):length(at) 
      Rij = pos[i] - pos[j]
      fqμ = force_q_μ(Rij, Qs[i], mus[j], λ)[1]
      fs[i] -= fqμ
      fs[j] += fqμ
      fμμ = force_μ_μ(Rij, mus[i], mus[j], λ)[1]
      fs[i] -= fμμ
      fs[j] += fμμ  
    end
    for j in Js
      Fs[j] += transpose(dV[j]) * fs[i]
    end
  end
  for (i, R) in enumerate(pos)
    for j = (i+1):length(Qs)
      Rij = R - pos[j]
      fqq = force_q_q(Rij,Qs[i], Qs[j], λ)[1]
      Fs[i] -= fqq
      Fs[j] += fqq
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
function electrostatic_energy(pos::AbstractArray, charges::AbstractVector, dipoles::AbstractArray, λ::Real, pbc::Bool=false)
  @assert pbc == false "Periodic boundary condition not yet supported"
  qq = 0
  qμ = 0
  μμ = 0
  for (i, R) in enumerate(pos)
    for j = (i+1):length(charges)
      Rij = R - pos[j]
      qq += soft_coulomb(Rij, charges[i], charges[j], λ)
      qμ += soft_q_μ(Rij, charges[i], dipoles[j], λ)
      μμ += soft_μ_μ(Rij, dipoles[i], dipoles[j], λ)
    end
  end
  return qq + qμ + μμ
end

"""
Electrostatic forces of a set of fixed point charges and fixed point dipoles
calculated using the soft core potentials. 
λ = 1.0 returns the non-soft core version.
"""
function electrostatic_forces(pos::AbstractArray, charges::AbstractVector, dipoles::AbstractArray, λ::Real, pbc::Bool=false)
  @assert pbc == false "Periodic boundary condition not yet supported"
  Fs = zeros(JVec{Float64}, length(charges))
  for (i, R) in enumerate(pos)
    for j = (i+1):length(charges)
      Rij = R - pos[j]
      fqq = force_q_q(Rij,charges[i], charges[j], λ)[1]
      Fs[i] -= fqq
      Fs[j] += fqq
      fqμ = force_q_μ(Rij, charges[i], dipoles[j], λ)[1]
      Fs[i] -= fqμ
      Fs[j] += fqμ
      fμμ = force_μ_μ(Rij, dipoles[i], dipoles[j], λ)[1]
      Fs[i] -= fμμ
      Fs[j] += fμμ      
    end
  end
  return Fs
end

"""
Soft core Coulomb interaction between two point charges
If λ=1.0 the normal Coulomb is returned
Adapted from LAMMPS `pair_style coul/long/soft`
"""

# TODO one single r_ij argument for forces
function soft_coulomb(Rij::AbstractVector, q1::Real, q2::Real, λ::Real=0.0, α::Real=10.0)
  r = norm(Rij)
  return q1 * q2 / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^0.5) * e * 1e10
end

function force_q_q(Rij::AbstractVector, q1::Real, q2::Real, λ::Real=0.0, α::Real=10.0)
  return Zygote.gradient(r -> soft_coulomb(r, q1, q2, λ, α), Rij)
end

"""
Soft core Coulomb interaction between point charge and point dipole
If λ=1.0 the normal Coulomb is returned
Analogously defined to LAMMPS `pair_style coul/long/soft`
"""
function soft_q_μ(Rij::AbstractVector, q1::Real, μ::AbstractArray, λ::Real=0.0, α::Real=10.0)
  #r = norm(Rij)
  #return q1 * sum(Rij .* μ) / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5) *1e10^2 * (1e-21/c_light)
  return q1 * dot(Rij, μ) / (4*π*ϵ_0 * (α * (1 - λ)^2 + dot(Rij, Rij))^1.5) *1e10^2 * (1e-21/c_light)
end

function force_q_μ(Rij::AbstractVector, q1::Real, μ::AbstractArray, λ::Real=0.0, α::Real=10.0)
  return Zygote.gradient((r, mu) -> soft_q_μ(r, q1, mu, λ, α), Rij, μ)
end

"""
Soft core Coulomb interaction between two point dipoles
If λ=1.0 the normal Coulomb is returned
Analogously defined to LAMMPS `pair_style coul/long/soft`
"""
function soft_μ_μ(Rij::AbstractVector, μ1::AbstractArray, μ2::AbstractArray, λ::Real=0.0, α::Real=10.0)
  #r = norm(Rij)
  #return 1 / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5) * (sum(μ1 .* μ2) - 3 * (sum(μ1 .* Rij)) * (sum(Rij .* μ2)) / (α * (1 - λ)^2 + r^2)) *1e10^3 * (1e-21/c_light)^2 / e
  return 1 / (4*π*ϵ_0 * (α * (1 - λ)^2 + dot(Rij, Rij))^1.5) * (dot(μ1, μ2) - 3 * dot(μ1, Rij) * dot(Rij, μ2) / (α * (1 - λ)^2 + dot(Rij, Rij))) *1e10^3 * (1e-21/c_light)^2 / e
end

function force_μ_μ(Rij::AbstractVector, μ1::AbstractArray, μ2::AbstractArray, λ::Real=0.0, α::Real=10.0)
  return Zygote.gradient((r, mu1, mu2) -> soft_μ_μ(r, mu1, mu2, λ, α), Rij, μ1, μ2)
end

end # end of module Electrostatics