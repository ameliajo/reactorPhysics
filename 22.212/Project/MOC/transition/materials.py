import numpy as np
from main import NGROUP

def import_xs(ngroup):
   """Create a MATERIALS dictionary using files output by an OpenMC script
   Parameters
   ----------
   ngroup : int
       Number of energy groups
   
   Returns
   -------
   MATERIALS : dict
       dictionary containing the cross section information for materials in the
       problem
   """
   folder = './make_xs/xs_'+str(ngroup)+'group/'
   fuel_file = '_cell_1'
   mod_file = '_cell_0'
   MATERIALS = {'fuel': {'total': np.loadtxt(folder+'total'+fuel_file),
                         'nufission': np.loadtxt(folder+'nufission'+fuel_file),
                         'scatter': np.loadtxt(folder+'scatter'+fuel_file),
                         'chi': np.loadtxt(folder+'chi'+fuel_file)
                         },
                'mod': {'total': np.loadtxt(folder+'total'+mod_file),
                        'nufission': np.loadtxt(folder+'nufission'+mod_file),
                        'scatter': np.loadtxt(folder+'scatter'+mod_file),
                        'chi': np.loadtxt(folder+'chi'+mod_file)
                       }
               }
   return MATERIALS

MATERIALS = import_xs(NGROUP)
