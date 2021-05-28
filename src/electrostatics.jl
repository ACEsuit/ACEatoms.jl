"""
`module Electro`

This module implements various functions related to (long range) electrostatic interactions
"""

module Electro

using LinearAlgebra

export get_dipole, electrostatic_energy

const μ_0 = 4.0e-7 * π
const c = 299792458.0
const ϵ_0 = 1 / μ_0 / c^2

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
  return sum(pos .* charges .+ dipoles, dims=1)

end

function electrostatic_energy(pos::Array, charges::Array, dipoles::Array, pbc::Bool=false)
  """
    Total electrostatic enenrgy of a set of point charges and point dipoles
    calculated using the soft core potentials. 
  """
  @assert pbc == false "Periodic boundary condition not yet supported"
  qq = 0
  qμ = 0
  μμ = 0
  for (i, R) in enumerate(pos)
    for j = (i+1):length(charges)
      qq += soft_coulomb(R, charges[i], pos[j], charges[j])
      qμ += soft_q_μ(R, charges[i], pos[j], dipoles[j])
      μμ += soft_μ_μ(R, dipoles[i], pos[j], dipoles[j])
    end
  end
  return qq + qμ + μμ
end

function soft_coulomb(pos1::Array, q1::Type, pos2::Array, q2::Type, λ::Type=0.9, α::Type=10.0)
  """
    Soft core Coulomb interaction between two point charges
    If λ=1.0 the normal Coulomb is returned
    Taken from LAMMPS `pair_style coul/long/soft`
  """
  r = norm(pos1 - pos2)
  return λ * q1 * q2 / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^0.5)
end

function soft_q_μ(pos1::Array, q1::Type, pos2::Array, μ::Array, λ::Type=0.9, α::Type=10.0)
  """
    Soft core Coulomb interaction between point charge and point dipole
    If λ=1.0 the normal Coulomb is returned
    Analogously defined to LAMMPS `pair_style coul/long/soft`
  """
  r12 = pos1 - pos2
  r = norm(r12)
  return λ * q1 * r12 ⋅ μ / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5)
end

function soft_μ_μ(pos1::Array, μ1::Type, pos2::Array, μ2::Array, λ::Type=0.9, α::Type=10.0)
  """
    Soft core Coulomb interaction between two point dipoles
    If λ=1.0 the normal Coulomb is returned
    Analogously defined to LAMMPS `pair_style coul/long/soft`
  """
  r12 = pos1 - pos2
  r = norm(r12)
  return λ / (4*π*ϵ_0 * (α * (1 - λ)^2 + r^2)^1.5) * (μ1 ⋅ μ2 - 3 * (μ1 ⋅ r12) * (r12 ⋅ μ2) / (α * (1 - λ)^2 + r^2))
end


end  # end of module Electro