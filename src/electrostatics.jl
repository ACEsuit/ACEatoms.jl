"""
`module Electrostatics`

This module implements various functions related to (long range) electrostatic interactions
"""

module Electrostatics

using LinearAlgebra
using Zygote
using JuLIP

export get_dipole, electrostatic_energy, soft_coulomb, 
        soft_q_μ, soft_μ_μ, soft_q_μ2, dipole, Fixed_q_dipole

const μ_0 = 4.0e-7 * π
const c = 299792458.0
const ϵ_0 = 1 / μ_0 / c^2
const e = 1.602176634e-19
const c_light = 299792458

"""
Get the overall dipole moment of a set of oint charges and point dipoles
Arguments:
  pos::Array natoms x 3 Array
  charges::Array natoms x 1 Array
  dipoles::Array natoms x 3 Array
Returns:
  total_dipole::Array 1 x 3 Array
"""
function get_dipole(pos::AbstractArray, charges::AbstractVector, dipoles::AbstractArray, pbc::Bool=false)
  @assert pbc == false "Periodic boundary condition not yet supported"
  return sum(pos .* charges .* (1e-11/c_light/e) .+ dipoles, dims=1)
end


#abstract type Fixed_q_dipole{T} end

struct Fixed_q_dipole{T} #<: Fixed_q_dipole{T} 
  x::T  
end

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
      qq += soft_coulomb(R, charges[i], pos[j], charges[j], λ)
      qμ += soft_q_μ(R, charges[i], pos[j], dipoles[j], λ)
      μμ += soft_μ_μ(R, dipoles[i], pos[j], dipoles[j], λ)
    end
  end
  return qq + qμ + μμ
end

"""
Soft core Coulomb interaction between two point charges
If λ=1.0 the normal Coulomb is returned
Adapted from LAMMPS `pair_style coul/long/soft`
"""
function soft_coulomb(pos1::AbstractVector, q1::Real, pos2::AbstractVector, q2::Real, λ::Real=0.9, α::Real=10.0)
  r = norm(pos1 - pos2)
  return λ * q1 * q2 / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^0.5) * e * 1e10
end

function force_q_q(pos1::AbstractVector, q1::Real, pos2::AbstractVector, q2::Real, λ::Real=0.9, α::Real=10.0)
  return gradient((p1, p2) -> soft_coulomb(p1, q1, p2, q2, λ, α), pos1, pos2)
end

"""
Soft core Coulomb interaction between point charge and point dipole
If λ=1.0 the normal Coulomb is returned
Analogously defined to LAMMPS `pair_style coul/long/soft`
"""
function soft_q_μ(pos1::AbstractVector, q1::Real, pos2::AbstractVector, μ::AbstractArray, λ::Real=0.9, α::Real=10.0)
  r12 = pos1 - pos2
  r = norm(r12)
  return λ * q1 * sum(r12 .* μ) / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5) *1e10^2 * (1e-21/c_light)
end

function force_q_μ(pos1::AbstractVector, q1::Real, pos2::AbstractVector, μ::AbstractArray, λ::Real=0.9, α::Real=10.0)
  return gradient((p1, p2) -> soft_q_μ(p1, q1, p2, μ, λ, α), pos1, pos2)
end

"""
Soft core Coulomb interaction between two point dipoles
If λ=1.0 the normal Coulomb is returned
Analogously defined to LAMMPS `pair_style coul/long/soft`
"""
function soft_μ_μ(pos1::AbstractVector, μ1::AbstractArray, pos2::AbstractVector, μ2::AbstractArray, λ::Real=0.9, α::Real=10.0)
  r12 = pos1 - pos2
  r = norm(r12)
  return λ / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5) * (sum(μ1 .* μ2) - 3 * (sum(μ1 .* r12)) * (sum(r12 .* μ2)) / (α * (1 - λ)^2 + r^2)) *1e10^3 * (1e-21/c_light)^2 / e 
end

function force_μ_μ(pos1::AbstractVector, μ1::AbstractArray, pos2::AbstractVector, μ2::AbstractArray, λ::Real=0.9, α::Real=10.0)
  return gradient((p1, p2) -> soft_μ_μ(p1, μ1, p2, μ2, λ, α), pos1, pos2)
end


end  # end of module Electrostatics