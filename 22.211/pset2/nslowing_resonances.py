import numpy as np
import matplotlib.pyplot as plt
from math import log
from res_construction import xs_from_res


def zerod_mc(mf_ratio,temp, neutrons, logemin, logemax, nbins, pick_res=[False,0], by_component=False):
    #Resonances used in this problem
    res_E = [6.673491e+0, 2.087152e+1, 3.668212e+1]
    J = [5.000000e-1, 5.000000e-1, 5.000000e-1]
    gn = [1.475792e-3, 1.009376e-2, 3.354568e-2]
    gg = [2.300000e-2, 2.286379e-2, 2.300225e-2]
    gfa = [0.000000, 5.420000e-8, 0.000000]
    gfb = [9.990000e-9, 0.000000, 9.770000e-9]
    if pick_res[0]:
        index = pick_res[1]
        res_E =  [res_E[index]]
        J     =  [J[index]]     
        gn    =  [gn[index]]    
        gg    =  [gg[index]]    
        gfa   =  [gfa[index]]   
        gfb   =  [gfb[index]]
    comp = 'all'
    ap = 0.948 #[barns]
    A = 238
    
    e_bins = np.logspace(0,2,nbins)
    e_freq = np.zeros(nbins-1)
    scat_freq = np.zeros(nbins-1)
    
    captures = 0
    for i in range(neutrons):
        energy = 1000 
        collisions = 0
        while energy > 1:
            e_freq_step = np.histogram(energy,e_bins)[0]
            mod_xs_micro = 20
            mod_xs_macro = mod_xs_micro*mf_ratio
            elastic_xs_micro = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,energy,reaction='elastic',comp=comp)
            elastic_xs_macro = elastic_xs_micro
            capture_xs_micro = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,energy,reaction='capture',comp=comp)
            capture_xs_macro = capture_xs_micro
            fuel_xs_macro = elastic_xs_macro+capture_xs_macro
            total_xs_macro = mod_xs_macro + fuel_xs_macro
            e_freq += e_freq_step/total_xs_macro
            mf_random = np.random.uniform(0,1)
            if mf_random < fuel_xs_macro/total_xs_macro:                         #Fuel reaction
                collision_random = np.random.uniform(0,1)
                if collision_random < capture_xs_macro/fuel_xs_macro:     #Capture in fuel
                    captures += 1
                    break
                scat_freq += e_freq_step
                alpha = ((A-1)/(A+1))**2
                energy = (energy-alpha*energy)*np.random.uniform(0,1)+alpha*energy    #Scatter in fuel
            elif mf_random > fuel_xs_macro/total_xs_macro:
                energy = np.random.uniform(0,1)*energy #Scatter in moderator
    e_freq = np.concatenate([[0],e_freq])/neutrons
    scat_freq = np.concatenate([[0],scat_freq])/neutrons
    abs_ratio = captures/neutrons
    if pick_res[0]:
        return e_bins,e_freq,abs_ratio,scat_freq
    elif by_component:
        elastic_xs_tot = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,reaction='elastic')
        capture_xs_tot = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,reaction='capture')
        xs_tot = elastic_xs_tot + capture_xs_tot

        elastic_xs_psi = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,reaction='elastic',comp='psi')
        psi_frac = elastic_xs_psi/xs_tot
        elastic_xs_chi = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,reaction='elastic',comp='chi')
        chi_frac = elastic_xs_chi/xs_tot
        elastic_xs_pot = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,reaction='elastic',comp='pot')
        pot_frac = elastic_xs_pot/xs_tot

        # Used to check the fractional components
        # plt.step(e_bins,psi_frac)
        # plt.step(e_bins,chi_frac)
        # plt.step(e_bins,pot_frac)
        # plt.step(e_bins,pot_frac+chi_frac+psi_frac)
        # plt.show()

        comp_wise = []
        comp_wise.append(e_freq*psi_frac)
        comp_wise.append(e_freq*chi_frac)
        comp_wise.append(e_freq*pot_frac)
        return e_bins, e_freq, abs_ratio, scat_freq, comp_wise
    else:
        return e_bins,e_freq,abs_ratio
