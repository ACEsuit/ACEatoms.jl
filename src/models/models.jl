module Models

import ACE 
ACE.@aceimports
ACE.@baseimports
ACE.@extimports

import ACEatoms: environment, AtomicEnvironment

import JuLIP

import JuLIP: SitePotential,
              cutoff,
              z2i, i2z, numz,
              AbstractCalculator,
              Atoms,
              chemical_symbol, AtomicNumber,
              JVec, JMat,
              energy, forces, virial

import JuLIP.MLIPs: IPBasis

import JuLIP.Potentials: ZList, SZList, zlist



include("twolayer.jl")

end