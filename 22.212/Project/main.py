from math import pi
import sys
sys.path.append('./dilutionTables/')

from interpret import *


class Nuclide:
    def __init__(self,N,r):
        self.N   = N
        self.pot = 4.0*pi*r*r
        self.sigT = []
        self.sigF = []
        self.sigA = []
        self.nuBar = []
        
def interpolate(x1,y1,x2,y2,x):
    return y1 + ((y2-y1)/(x2-x1))*x - (y2*x1-y1*x1)/(x2-x1)

def getDataFromEqTable(boundingDilutions,nuclideData):
    d1     = boundingDilutions[0][1]
    sigT1  = nuclideData.sigT[boundingDilutions[0][0]]
    sigF1  = nuclideData.sigF[boundingDilutions[0][0]]
    sigA1  = nuclideData.sigA[boundingDilutions[0][0]]
    nuBar1 = nuclideData.nuBar
    
    d2     = boundingDilutions[1][1]
    sigT2  = nuclideData.sigT[boundingDilutions[1][0]]
    sigF2  = nuclideData.sigF[boundingDilutions[1][0]]
    sigA2  = nuclideData.sigA[boundingDilutions[1][0]]
    nuBar2 = nuclideData.nuBar

    sigT = (interpolate(d1,sigT1,d2,sigT2,backgroundXS))
    sigF = (interpolate(d1,sigF1,d2,sigF2,backgroundXS))
    sigA = (interpolate(d1,sigA1,d2,sigA2,backgroundXS))
    nuBar = (interpolate(d1,nuBar1,d2,nuBar2,backgroundXS))

    return sigT,sigF,sigA,nuBar



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

# background = sum_nonResNuclides [N_nonRes * sig_pot_nonRes / N_res] + 1/(N_res*l_bar)
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


for g in range(nGroups):
    sigT,sigF,sigA,nuBar = getDataFromEqTable(boundingDilutions,groupsU235[g])
    u235.sigT.append(sigT)
    u235.sigF.append(sigF)
    u235.sigA.append(sigA)
    u235.nuBar.append(nuBar)

print(u235.sigT)
print(u235.sigF)
print(u235.sigA)
print(u235.nuBar)






###############################################################################
# Evaluate group-wise collision probability using the e ective cross section evaluated in 
###############################################################################

###############################################################################
# Update the background cross section 
###############################################################################




##############################################################################
# Repeat
##############################################################################


