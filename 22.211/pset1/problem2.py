import numpy as np
import matplotlib.pyplot as plt

energies = [ 6.673491E+0, 2.087152E+1, 3.668212E+1 ]
gamma_n  = [ 1.475792E-3, 1.009376E-2, 3.354568E-2 ]
gamma_g  = [ 2.300000E-2, 2.286379E-2, 2.300225E-2 ]
gamma    = [gamma_n[i]+gamma_g[i] for i in range(len(gamma_n))]
sigma_pot = 11.2934
A = 238.0

E_vec = np.linspace(1e-2,1e2,10000)
xs_gamma = [0.0]*len(E_vec)
xs_elast = [0.0]*len(E_vec)
xs_total = [0.0]*len(E_vec)

for i in range(len(E_vec)): 
    E = E_vec[i]
    for j in range(len(energies)):
        E0 = energies[j]
        x = 2.0*(E-E0)/gamma[j]
        psi = 1.0/(1.0+x*x)
        chi =  x /(1.0+x*x)

        r = 2603911/E0 * (A+1)/A
        q = 2.0*(r*sigma_pot)**0.5 

        xs_gamma[i] += gamma_n[j]*gamma_g[j]/(gamma[j]**2) * (E0/E)**0.5 * (r*psi)
        xs_elast[i] += gamma_n[j]*gamma_n[j]/(gamma[j]**2) * ( (r*psi) + gamma[j]/gamma_n[j]*q*chi) + sigma_pot
    xs_total[i] = xs_gamma[i] + xs_elast[i]

plt.plot(E_vec,xs_gamma)
plt.plot(E_vec,xs_elast)
plt.plot(E_vec,xs_total)
plt.xscale('log')
plt.yscale('log')
plt.show()





