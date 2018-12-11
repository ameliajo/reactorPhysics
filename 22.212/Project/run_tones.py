from math import pi
import matplotlib.pyplot as plt
from numpy import ma

import sys
sys.path.append('./dilutionTables/')
sys.path.append('./GettingXS/')
sys.path.append('./collisionProb/')

from interpret import *
from import_XS_help import *
from XS_nuclideSpecific import *
from getCollisionProb import *
from XS import modTotal0
import matplotlib.colors as colors
import matplotlib.cm as cmx



class Pin:
    def __init__(self,N,radius,all_nuclides_pin):
        self.U235 = Nuclide(N['U235'],radius['U235'],all_nuclides_pin[0])
        self.U238 = Nuclide(N['U238'],radius['U238'],all_nuclides_pin[1])
        self.O16  = Nuclide(N['O16' ],radius['O16' ],all_nuclides_pin[2])

    def calcSigT(self,nGroups):
        self.SigT = []
        for g in range(nGroups):
            self.SigT.append(self.U235.SigT[g]        + \
                             self.U238.openMC_SigT[g] + \
                             self.O16.openMC_SigT[g])




#jet = cm = plt.get_cmap('hot')
#jet = cm = plt.get_cmap('tab10')
#jet = cm = plt.get_cmap('autumn')
#cNorm  = colors.Normalize(vmin=0, vmax=8)
#scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# Define problem geometry
pinRad = 0.39128;  pitch = 1.26
l_bar  = 2*pinRad; C     = 0.15


E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]




###############################################################################
# Load in dilution table
###############################################################################

# groupsU235[g].sigT[d] will give you sigT for group g for dilution index d
groupsU235 = interpretGENDF('u235',False,'dilutionTables/u235_tape26')
groupsU235.reverse() # Because NJOY does groups in inverse order (low -> high)
dilution = groupsU235[0].dilutionVals # All dilution values are the same, it
                                      # doesn't matter I'm pulling this from 
                                      # energy group 0

# Gotten from OpenMC input, using uo2.get_nuclide_atom_densities and multiplying
# by 1E24 bc the number densities are provided in atoms / b cm
N_dict_Hi = { 'U235' : N_hi_U235, 'U238' : N_hi_U238, 'O16'  : N_hi_O16 }
N_dict_Lo = { 'U235' : N_lo_U235, 'U238' : N_lo_U238, 'O16'  : N_lo_O16 }

# Called ``AP'' in ENDF, found in MF 2 MT 151
radius = { 'U234': 8.930000E-1, 'U235': 9.602000E-1, \
           'U238': 9.480000E-1, 'U236': 9.354000E-1, \
           'O16' : 5.562563E-1, 'O17' : 5.780000E-1 }


# Define pins materials
pins = [Pin([N_dict_Hi,N_dict_Lo][i%2],radius,all_nuclides[i]) for i in range(9)]


###############################################################################
# Assume initial background cross sections for resonance nuclides.  
# They can be evaluated using the conventional equivalence method with the 
# Dancoff  correction.
###############################################################################

#------------------------------------------------------------------------------
# We're going to treat U-235 as the resonant nuclide for now, for pin0
#------------------------------------------------------------------------------
nonRes = [pins[0].U238,pins[0].O16]

# background = sum_nonRes [N_nonRes * sig_potNonRes / N_res] + 1/(N_res*l_bar)
sig0 = sum([nuclide.N*nuclide.pot for nuclide in nonRes])/pins[0].U235.N + \
       1.0 / (pins[0].U235.N*l_bar*(1.0-C))


nGroups = len(groupsU235)
sig0EnergyVec = [sig0]*nGroups

converged = False
counter = 0
sig0Vals = [sig0]
while not converged:
    print(counter)
    #print(sig0EnergyVec)

    ###########################################################################
    # Evaluate the effective cross sections of resonance nuclides using the 
    # conventional equivalence theory.
    ###########################################################################

    boundingDilutionsVec = [None]*nGroups
    for g in range(nGroups):
        for i in range(len(dilution)-1):
            if dilution[i] > sig0EnergyVec[g] > dilution[i+1]:
                boundingDilutionsVec[g] = [(i,dilution[i]),(i+1,dilution[i+1])]
                break

    if None in boundingDilutionsVec: raise ValueError('Ideal dilution value out of range')

    # Pulling cross sections from the dilution table
    for pin in pins:
        pin.U235.resetXS()
        for g in range(nGroups):
            pin.U235.addXS(getDataFromEqTable(boundingDilutionsVec[g],\
                           groupsU235[g],sig0EnergyVec[g]))

        # Making these microscopic cross sections into macroscropic cross sections
        pin.U235.convertToMacro()



    ###########################################################################
    # Evaluate group-wise collision probability using the effctive XS from dilution
    ###########################################################################

    #--------------------------------------------------------------------------
    # Calculate the macro. SigT for the high/low enr. fuel pins, and plug into the
    # collision probability MC script
    #--------------------------------------------------------------------------
    SigT_hi = [ sum([nucl.openMC_SigT[g] for nucl in nonRes]) + \
                pins[0].U235.SigT[g] for g in range(nGroups) ]
    SigT_lo = [ sum([nucl.openMC_SigT[g] for nucl in nonRes]) + \
                pins[1].U235.SigT[g] for g in range(nGroups) ]

    # For hi/lo enr. of fuel, use values we just pulled from dilution table
    # For moderator, use openMC values generated from grid_3x3.py
    collisionProbs = [                                                   \
        getCollisionProb( pitch, pinRad, plot=False, numParticles=2000,  \
          hole=False, fSigT_hi=SigT_hi[g], fSigT_lo=SigT_lo[g],          \
          mSigT=modTotal0[g], verbose=False, startNeutronsFrom=0 )       \
        for g in range(nGroups)]





    ###########################################################################
    # Update the background cross section 
    ###########################################################################

    # sig0 = SUM_pins SUM_nonRes P(pin->me) * V(pin) * N(nonResInPin)*sig(pot)
    #        ------------------------------------------------------------------
    #                 SUM_pins * P(pin->me) * V(pin) * N(resInPin)

    tones_Numer = 0.0
    tones_Denom = 0.0

    for pinID,pin in enumerate(pins):
        pin.calcSigT(nGroups)
        P_0_to_i = np.array([collisionProbs[g][pinID] for g in range(nGroups)])
        # RECIPROCITY RELATION
        # P(i->0) = P(0->i) * SigT_0 / SigT_i
        P_i_to_0 = P_0_to_i * pins[0].SigT / pin.SigT

        SUM_nonRes_sigPot = pin.U238.N*pin.U238.pot*1E-24 +\
                            pin.O16.N*pin.O16.pot*1E-24
        
        tones_Numer += P_i_to_0 * SUM_nonRes_sigPot
        tones_Denom += P_i_to_0 * pin.U235.N

    sig0_new = sum(tones_Numer/tones_Denom*1e24)/nGroups
    #plt.step(E_bounds,[1e24*tones_Numer[0]/tones_Denom[0]]+list(1e24*tones_Numer/tones_Denom),label='inter'+str(counter),color=scalarMap.to_rgba(counter))


    ###########################################################################
    # Repeat
    ###########################################################################
    sig0 = sig0_new 
    sig0EnergyVec = tones_Numer/tones_Denom*1e24

    counter += 1
    if counter > 5:
        print(counter)
        SigT_hi = [ sum([nucl.openMC_SigT[g] for nucl in nonRes]) + \
                    pins[0].U235.SigT[g] for g in range(nGroups) ]
        SigT_lo = [ sum([nucl.openMC_SigT[g] for nucl in nonRes]) + \
                    pins[1].U235.SigT[g] for g in range(nGroups) ]
        SigA_hi = [ sum([nucl.openMC_SigA[g] for nucl in nonRes]) + \
                    pins[0].U235.SigA[g] for g in range(nGroups) ]
        SigA_lo = [ sum([nucl.openMC_SigA[g] for nucl in nonRes]) + \
                    pins[1].U235.SigA[g] for g in range(nGroups) ]
        nuSigF_hi = [ sum([nucl.openMC_nuSigF[g] for nucl in nonRes]) + \
                    pins[0].U235.SigF[g]*pins[0].U235.nuBar[g] for g in range(nGroups) ]
        nuSigF_lo = [ sum([nucl.openMC_nuSigF[g] for nucl in nonRes]) + \
                    pins[1].U235.SigF[g]*pins[0].U235.nuBar[g] for g in range(nGroups) ]

        f = open("XS-tones.py","a")
        for i in range(9):
            f.write("fuelTotal"+str(i)+"  = "+str([float("%.8f"%[SigT_hi,SigT_lo][i%2][g]) for g in range(nGroups)])+"\n")
            f.write("fuelAbsorption"+str(i)+"  = "+str([float("%.8f"%[SigA_hi,SigA_lo][i%2][g]) for g in range(nGroups)])+"\n")
            f.write("fuelNuFission"+str(i)+"  = "+str([float("%.8f"%[nuSigF_hi,nuSigF_lo][i%2][g]) for g in range(nGroups)])+"\n")

            f.write("\n\n")
        f.close()

        frac = []
        for i in range(10):
            print(sig0EnergyVec[i])
            print(pins[0].U235.pot)
            print(pins[0].U235.pot)
            numer = sig0EnergyVec[i]+pins[0].U235.pot
            denom = sig0EnergyVec[i]+pins[0].U235.sigT[i]
            frac.append(numer/denom)

        print(frac)

        break

"""
ax = plt.gca()
ax.set_facecolor('xkcd:pale grey')
ax.set_facecolor('xkcd:off white')
plt.xscale('log')
plt.xlabel('Energy (eV)')
plt.ylabel('Sig0 approximation')
plt.legend(loc='best')
#plt.savefig('sig0Estimations.png')
plt.show()

"""



