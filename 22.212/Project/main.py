from math import pi
import sys
sys.path.append('./dilutionTables/')
sys.path.append('./GettingXS/')

from interpret import *
from import_XS_help import *
from XS_nuclideSpecific import *




# Load in dilution table
groupsU238 = interpretGENDF('u238',False,'dilutionTables/u238_tape26')
groupsU235 = interpretGENDF('u235',False,'dilutionTables/u235_tape26')

# Gotten from OpenMC input, using uo2.get_nuclide_atom_densities and multiplying
# by 1E24 bc the number densities are provided in atoms / b cm
numberDensities = { 'U234': 6.68004682609e+18, 'U235': 7.47364244574e+20, \
                    'U238': 2.23122440162e+22, 'U236': 3.42328721032e+18, \
                    'O16' : 4.61219363482e+22,  'O17' : 1.74868413889e+19 }

# Called ``AP'' in ENDF, found in MF 2 MT 151
radius = { 'U234': 8.930000E-1, 'U235': 9.602000E-1, \
           'U238': 9.480000E-1, 'U236': 9.354000E-1, \
           'O16' : 5.562563E-1, 'O17' : 5.780000E-1 }


u234 = Nuclide(numberDensities['U234'],radius['U234'])
u235 = Nuclide(numberDensities['U235'],radius['U235'])
u236 = Nuclide(numberDensities['U236'],radius['U236'])
u238 = Nuclide(numberDensities['U238'],radius['U238'])
o16  = Nuclide(numberDensities['O16' ],radius['O16' ])
o17  = Nuclide(numberDensities['O17' ],radius['O17' ])
u234.importFromOpenMC(u234_all)
u235.importFromOpenMC(u235_all)
u236.importFromOpenMC(u236_all)
u238.importFromOpenMC(u238_all)
o16.importFromOpenMC(o16_all)
o17.importFromOpenMC(o17_all)


l_bar = 2.0 * 0.39128
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
backgroundXS = 0.0
for nuclide in nonRes:
    backgroundXS += nuclide.N*nuclide.pot
backgroundXS /= res.N
backgroundXS += 1.0 / (res.N*l_bar*(1.0-C))


###############################################################################
# Evaluate the effective cross sections of resonance nuclides using the 
# conventional equivalence theory.
###############################################################################

dilution = groupsU235[0].dilutionVals
nGroups = len(groupsU235)
boundingDilutions = None
for i in range(len(dilution)-1):
    if dilution[i] > backgroundXS > dilution[i+1]:
        boundingDilutions = [(i,dilution[i]),(i+1,dilution[i+1])]
        break

if boundingDilutions == None: raise ValueError('Ideal dilution value out of range')


# Pulling cross sections from the dilution table
for g in range(nGroups):
    u235.addXS(getDataFromEqTable(boundingDilutions,groupsU235[g],backgroundXS))


# Making these microscopic cross sections into macroscropic cross sections
u235.convertToMacro()



###############################################################################
# Evaluate group-wise collision probability using the effctive XS from dilution
###############################################################################
SigT_new = [0.0]*nGroups
for nuclide in nonRes:
    SigT_new = [SigT_new[g] + nuclide.openMC_SigT[g] for g in range(nGroups)]
SigT_new = [SigT_new[g] + u235.SigT[g] for g in range(nGroups)]
    
print(SigT_new)









###############################################################################
# Update the background cross section 
###############################################################################




##############################################################################
# Repeat
##############################################################################


