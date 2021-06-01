"""
`module Electrostatics`

This module implements various functions related to (long range) electrostatic interactions
"""

module Electrostatics

using LinearAlgebra

export get_dipole, electrostatic_energy, soft_coulomb, soft_q_μ, soft_μ_μ

const μ_0 = 4.0e-7 * π
const c = 299792458.0
const ϵ_0 = 1 / μ_0 / c^2
const e = 1.602176634e-19
const c_light = 299792458


function get_dipole(pos::Array, charges::Array, dipoles::Array, pbc::Bool=false)
  """
  Get the overall dipole moment of a set of oint charges and point dipoles
  Arguments:
    pos::Array natoms x 3 Array
    charges::Array natoms x 1 Array
    dipoles::Array natoms x 3 Array
  Returns:
    total_dipole::Array 1 x 3 Array
  """
  @assert pbc == false "Periodic boundary condition not yet supported"
  return sum(pos .* charges .* (1e-11/c_light/e) .+ dipoles, dims=1)
end

function electrostatic_energy(pos::Array, charges::Array, dipoles::Array, λ::Any, pbc::Bool=false)
  """
    Total electrostatic enenrgy of a set of point charges and point dipoles
    calculated using the soft core potentials. 
    λ = 1.0 returns the non-soft core version.
  """
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

function soft_coulomb(pos1, q1, pos2, q2, λ=0.9, α=10.0)
  """
    Soft core Coulomb interaction between two point charges
    If λ=1.0 the normal Coulomb is returned
    Adapted from LAMMPS `pair_style coul/long/soft`
  """
  r = norm(pos1 - pos2)
  return λ * q1 * q2 / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^0.5) * e * 1e10
end

function soft_q_μ(pos1, q1, pos2, μ, λ=0.9, α=10.0)
  """
    Soft core Coulomb interaction between point charge and point dipole
    If λ=1.0 the normal Coulomb is returned
    Analogously defined to LAMMPS `pair_style coul/long/soft`
  """
  r12 = pos1 - pos2
  r = norm(r12)
  return λ * q1 * r12 ⋅ μ / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5) *1e10^2 * (1e-21/c_light)
end

function soft_μ_μ(pos1, μ1, pos2, μ2, λ=0.9, α=10.0)
  """
    Soft core Coulomb interaction between two point dipoles
    If λ=1.0 the normal Coulomb is returned
    Analogously defined to LAMMPS `pair_style coul/long/soft`
  """
  r12 = pos1 - pos2
  r = norm(r12)
  return λ / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5) * (μ1 ⋅ μ2 - 3 * (μ1 ⋅ r12) * (r12 ⋅ μ2) / (α * (1 - λ)^2 + r^2)) *1e10^3 * (1e-21/c_light)^2 / e 
end


end  # end of module Electro