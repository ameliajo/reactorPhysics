from math import pi
import sys
sys.path.append('./dilutionTables/')
sys.path.append('./GettingXS/')
sys.path.append('./collisionProb/')

from interpret import *
from import_XS_help import *
from XS_nuclideSpecific import *
from getCollisionProb import *
from XS import modTotal0


class Pin:
    def __init__(self,N,radius,all_nuclides_pin):
        self.nuclides = { \
        'U235' : Nuclide(N['U235'],radius['U235'],all_nuclides_pin[0]), \
        'U238' : Nuclide(N['U238'],radius['U238'],all_nuclides_pin[1]), \
        'O16'  : Nuclide(N['O16' ],radius['O16' ],all_nuclides_pin[2]) }

    def calcSigT(self,nGroups):
        self.SigT = []
        for g in range(nGroups):
            self.SigT.append(self.nuclides['U235'].SigT[g]        + \
                             self.nuclides['U238'].openMC_SigT[g] + \
                             self.nuclides['O16'].openMC_SigT[g])





# Load in dilution table
# groupsU235[g].sigT[d] will give you sigT for group g for dilution index d
groupsU235 = interpretGENDF('u235',False,'dilutionTables/u235_tape26')
groupsU235.reverse() # Because NJOY does groups in inverse order (low -> high)

# Gotten from OpenMC input, using uo2.get_nuclide_atom_densities and multiplying
# by 1E24 bc the number densities are provided in atoms / b cm
N_dict_Hi = { 'U235' : N_hi_U235, 'U238' : N_hi_U238, 'O16'  : N_hi_O16 }
N_dict_Lo = { 'U235' : N_lo_U235, 'U238' : N_lo_U238, 'O16'  : N_lo_O16 }

# Called ``AP'' in ENDF, found in MF 2 MT 151
radius = { 'U234': 8.930000E-1, 'U235': 9.602000E-1, \
           'U238': 9.480000E-1, 'U236': 9.354000E-1, \
           'O16' : 5.562563E-1, 'O17' : 5.780000E-1 }

pin_0 = Pin(N_dict_Hi,radius,all_nuclides_pin_0)
pin_2 = Pin(N_dict_Hi,radius,all_nuclides_pin_2)
pin_4 = Pin(N_dict_Hi,radius,all_nuclides_pin_4)
pin_6 = Pin(N_dict_Hi,radius,all_nuclides_pin_6)
pin_8 = Pin(N_dict_Hi,radius,all_nuclides_pin_8)
highPins = [pin_0,pin_2,pin_4,pin_6,pin_8]

pin_1 = Pin(N_dict_Lo,radius,all_nuclides_pin_1)
pin_3 = Pin(N_dict_Lo,radius,all_nuclides_pin_3)
pin_5 = Pin(N_dict_Lo,radius,all_nuclides_pin_5)
pin_7 = Pin(N_dict_Lo,radius,all_nuclides_pin_7)
lowPins = [pin_1,pin_3,pin_5,pin_7]
pins = highPins + lowPins


u235 = pin_0.nuclides['U235']
u238 = pin_0.nuclides['U238']
o16  = pin_0.nuclides['O16']

pinRadius = 0.39128
pitch = 1.26

l_bar = 2.0 * pinRadius
C = 0.15



###############################################################################
# Assume initial background cross sections for resonance nuclides.  
# They can be evaluated using the conventional equivalence method with the 
# Dancoff  correction.
###############################################################################

#------------------------------------------------------------------------------
# We're going to treat U-235 as the resonant nuclide for now, for pin0
#------------------------------------------------------------------------------
res    = u235
nonRes = [u238,o16]

# background = sum_nonRes [N_nonRes * sig_potNonRes / N_res] + 1/(N_res*l_bar)
sig0 = sum([nuclide.N*nuclide.pot for nuclide in nonRes])/res.N + 1.0/(res.N*l_bar*(1.0-C))


###############################################################################
# Evaluate the effective cross sections of resonance nuclides using the 
# conventional equivalence theory.
###############################################################################

dilution = groupsU235[0].dilutionVals
nGroups = len(groupsU235)
boundingDilutions = None
for i in range(len(dilution)-1):
    if dilution[i] > sig0 > dilution[i+1]:
        boundingDilutions = [(i,dilution[i]),(i+1,dilution[i+1])]
        break

if boundingDilutions == None: raise ValueError('Ideal dilution value out of range')


# Pulling cross sections from the dilution table
for pin in pins:
    for g in range(nGroups):
        pin.nuclides['U235'].addXS(getDataFromEqTable(boundingDilutions,groupsU235[g],sig0))

    # Making these microscopic cross sections into macroscropic cross sections
    pin.nuclides['U235'].convertToMacro()



###############################################################################
# Evaluate group-wise collision probability using the effctive XS from dilution
###############################################################################
SigT_new_hi = [0.0]*nGroups
SigT_new_lo = [0.0]*nGroups
for nuclide in nonRes:
    SigT_new_hi = [SigT_new_hi[g] + nuclide.openMC_SigT[g] for g in range(nGroups)]
    SigT_new_lo = [SigT_new_lo[g] + nuclide.openMC_SigT[g] for g in range(nGroups)]

SigT_new_hi = [SigT_new_hi[g] + pin_0.nuclides['U235'].SigT[g] for g in range(nGroups)]
SigT_new_lo = [SigT_new_lo[g] + pin_1.nuclides['U235'].SigT[g] for g in range(nGroups)]


nuclides = nonRes + [res]

collisionProbsFromCorner = []
collisionProbsFromSide   = []
collisionProbsFromCenter = []

for g in range(nGroups):
    fSigT_hi = SigT_new_hi[g]  # Use the values we just pulled from dilution table
    fSigT_lo = SigT_new_lo[g]  
    mSigT = modTotal0[g] # Use openMC values generated from grid_3x3.py

    collProb_from_corner = getCollisionProb( pitch, pinRadius, plot=False,   \
      numParticles=500, hole=False, fSigT_hi=fSigT_hi, fSigT_lo=fSigT_lo,    \
      mSigT=mSigT, verbose=False, startNeutronsFrom=0)

    collisionProbsFromCorner.append(collProb_from_corner)


#print(collisionProbsFromCorner[0])
#print(collisionProbsFromCorner[1])
#print(collisionProbsFromCorner[2])


P_0_to_0 = np.array([collisionProbsFromCorner[g][0] for g in range(nGroups)])
P_0_to_1 = np.array([collisionProbsFromCorner[g][1] for g in range(nGroups)])

#print(P_0_to_0)
print(P_0_to_1)
print()
pin_0.calcSigT(nGroups)
pin_1.calcSigT(nGroups)
pin_2.calcSigT(nGroups)
pin_3.calcSigT(nGroups)
pin_4.calcSigT(nGroups)
pin_5.calcSigT(nGroups)
pin_6.calcSigT(nGroups)
pin_7.calcSigT(nGroups)
pin_8.calcSigT(nGroups)

# P(1->0) = P(0->1) * SigT_0 / SigT_1
P_1_to_0 = P_0_to_1 * pin_0.SigT / pin_1.SigT
print(P_1_to_0)






###############################################################################
# Update the background cross section 
###############################################################################

# sig0 = SUM_pins SUM_nonRes P(pin->myPin) * V(pin) * N(nonRes in pin)*sig(pot)
#        ----------------------------------------------------------------------
#                 SUM_pins * P(pin->myPin) * V(pin) * N(res in pin)









##############################################################################
# Repeat
##############################################################################


