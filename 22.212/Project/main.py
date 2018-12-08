from math import pi
import sys
sys.path.append('./dilutionTables/')

from interpret import *


       

# Load in dilution table
groupsU238 = interpretGENDF('u238',False,'dilutionTables/u238_tape26')
groupsU235 = interpretGENDF('u235',False,'dilutionTables/u235_tape26')

# N = rho*N_A/M
N_A = 6.022E23
rho = 10.341 # g/cm^3

# Gotten from OpenMC input, using uo2.get_nuclide_atom_densities and multiplying
# by 1E24 bc the number densities are provided in atoms / b cm
numberDensities = { 'U234': 6.68004682609e+18, 'U235': 7.47364244574e+20, \
                    'U238': 2.23122440162e+22, 'U236': 3.42328721032e+18, \
                    'O16' : 4.61219363482e+22,  'O17' : 1.74868413889e+19 }

# Called ``AP'' in ENDF, found in MF 2 MT 151
radius = { 'U234': 8.930000E-1, 'U235': 9.602000E-1, \
           'U238': 9.480000E-1, 'U236': 9.354000E-1, \
           'O16' : 5.562563E-1, 'O17' : 5.780000E-1 }

potentialXS = {}
for nuclide in radius:
    potentialXS[nuclide] = 4.0*pi*radius[nuclide]**2

l_bar = 2.0 * 0.39128

# We're going to treat U-235 as the resonant nuclide for now

# background = sum_nonResNuclides [N_nonRes * sig_pot_nonRes / N_res] + 1/(N_res*l_bar)
backgroundXS = 0.0
for 

# Assume initial background cross sections for resonance nuclides.  They can be evaluated using the conventional equivalence method with the Danco  correction.

# Evaluate the e ective cross sections of resonance nuclides using the conventional equiva- lence theory.
# Evaluate group-wise collision probability using the e ective cross section evaluated in 
# Update the background cross section 
# Repeat


