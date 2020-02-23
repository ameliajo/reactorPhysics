import numpy as np
from math import pi
from scipy.special import wofz

def xs_from_res(ap, A, res_E, gn, gg, gfa, gfb, temp, energy, reaction=None, comp='all'):
    sigma_pot = 4*pi*ap**2
    k_b = 8.6173e-5 #[ev/K]
    xs = 0
    for i in range(len(res_E)):
        g_tot = gn[i]+gg[i]+gfa[i]+gfb[i]
        r = 2603911/res_E[i]*(A+1)/A
        q = 2*(r*sigma_pot)**0.5
        x = 2*(energy-res_E[i])/g_tot
        if temp==0:
            psi = 1/(1+x**2)
            chi = x/(1+x**2)
        else:
            xi = g_tot*(A/(4*k_b*temp*res_E[i]))**0.5
            fadeeva_eval = xi*wofz((x+1j)*xi)
            psi = pi**0.5*np.real(fadeeva_eval)
            chi = pi**0.5*np.imag(fadeeva_eval)
            if comp=='psi':
                chi = 0
            if comp=='chi':
                psi = 0
            if comp=='pot':
                chi = 0
                psi = 0
        if reaction=='capture':
            res_xs = gn[i]*gg[i]/g_tot**2*(res_E[i]/energy)**0.5*r*psi
        elif reaction=='elastic':
            res_xs = (gn[i]/g_tot)**2*(r*psi+g_tot/gn[i]*q*chi)
        xs += res_xs
    if reaction=='elastic' and comp!='psi' and comp!='chi':
        xs += sigma_pot
    return xs
