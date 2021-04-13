
# -----------------------------------------------------------------------------
# JuLIP, ACE, etc : just use import throughout, this avoids bugs


import Base: ==, length

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


import ACEbase: alloc_temp, alloc_temp_d,
                evaluate, evaluate_d,
                evaluate!, evaluate_d!,
                read_dict, write_dict,
                fltype, rfltype,
                alloc_B, alloc_dB,
                combine,
                ScalarACEBasis, _allfieldsequal

import ACE: AbstractState,
            AbstractDiscreteState,
            AbstractContinuousState,
            AbstractConfiguration,
            scaling 
