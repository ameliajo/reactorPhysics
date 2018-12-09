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




# Load in dilution table
# groupsU235[g].sigT[d] will give you sigT for group g for dilution index d
groupsU238 = interpretGENDF('u238',False,'dilutionTables/u238_tape26')
groupsU235 = interpretGENDF('u235',False,'dilutionTables/u235_tape26')
groupsU235.reverse() # Because NJOY does groups in inverse order (low -> high)

# Gotten from OpenMC input, using uo2.get_nuclide_atom_densities and multiplying
# by 1E24 bc the number densities are provided in atoms / b cm
numberDensities = { 'U234': 6.68004682609e+18, 'U235': 7.47364244574e+20, \
                    'U238': 2.23122440162e+22, 'U236': 3.42328721032e+18, \
                    'O16' : 4.61219363482e+22,  'O17' : 1.74868413889e+19 }

# Called ``AP'' in ENDF, found in MF 2 MT 151
radius = { 'U234': 8.930000E-1, 'U235': 9.602000E-1, \
           'U238': 9.480000E-1, 'U236': 9.354000E-1, \
           'O16' : 5.562563E-1, 'O17' : 5.780000E-1 }


u234 = Nuclide(numberDensities['U234'],radius['U234'],u234_all)
u235 = Nuclide(numberDensities['U235'],radius['U235'],u235_all)
u236 = Nuclide(numberDensities['U236'],radius['U236'],u236_all)
u238 = Nuclide(numberDensities['U238'],radius['U238'],u238_all)
o16  = Nuclide(numberDensities['O16' ],radius['O16' ],o16_all)
o17  = Nuclide(numberDensities['O17' ],radius['O17' ],o17_all)

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
# We're going to treat U-235 as the resonant nuclide for now
#------------------------------------------------------------------------------
res    = u235
nonRes = [u234,u236,u238,o16,o17]

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
for g in range(nGroups):
    u235.addXS(getDataFromEqTable(boundingDilutions,groupsU235[g],sig0))


# Making these microscopic cross sections into macroscropic cross sections
u235.convertToMacro()



###############################################################################
# Evaluate group-wise collision probability using the effctive XS from dilution
###############################################################################
SigT_new = [0.0]*nGroups
for nuclide in nonRes:
    SigT_new = [SigT_new[g] + nuclide.openMC_SigT[g] for g in range(nGroups)]
SigT_new = [SigT_new[g] + u235.SigT[g] for g in range(nGroups)]


nuclides = nonRes + [res]

collisionProbsFromCorner = []
collisionProbsFromSide   = []
collisionProbsFromCenter = []

for g in range(nGroups):
    fSigT = SigT_new[g]  # Use the values we just pulled from dilution table
    mSigT = modTotal0[g] # Use openMC values generated from grid_3x3.py

    collProb_from_corner = getCollisionProb( pitch, pinRadius, plot=False,   \
      numParticles=500, hole=False, fSigT=fSigT, mSigT=mSigT, verbose=False, \
      startNeutronsFrom=0)
    collProb_from_side = getCollisionProb( pitch, pinRadius, plot=False,     \
      numParticles=500, hole=False, fSigT=fSigT, mSigT=mSigT, verbose=False, \
      startNeutronsFrom=1)
    collProb_from_center = getCollisionProb( pitch, pinRadius, plot=False,   \
      numParticles=500, hole=False, fSigT=fSigT, mSigT=mSigT, verbose=False, \
      startNeutronsFrom=4)

    collisionProbsFromCorner.append(collProb_from_corner)
    collisionProbsFromSide.append(collProb_from_side)
    collisionProbsFromCenter.append(collProb_from_center)


print(collisionProbsFromCorner[0])
print(collisionProbsFromSide[0])
print(collisionProbsFromCenter[0])









###############################################################################
# Update the background cross section 
###############################################################################

# sig0 = SUM_pins SUM_nonRes P(pin->myPin) * V(pin) * N(nonRes in pin)*sig(pot)
#        ----------------------------------------------------------------------
#                 SUM_pins * P(pin->myPin) * V(pin) * N(res in pin)


P







##############################################################################
# Repeat
##############################################################################


