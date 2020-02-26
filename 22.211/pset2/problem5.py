import numpy as np
import matplotlib.pyplot as plt
from math import log
from res_construction import xs_from_res

def zerod_mc(mf_ratio,temp, neutrons, logemin, logemax, nbins, pick_res=[False,0], by_component=False):
    #Resonances used in this problem
    res_E = [6.673491e+0, 2.087152e+1, 3.668212e+1]
    gn = [1.475792e-3, 1.009376e-2, 3.354568e-2]
    gg = [2.300000e-2, 2.286379e-2, 2.300225e-2]
    gfa = [0.000000, 5.420000e-8, 0.000000]
    gfb = [9.990000e-9, 0.000000, 9.770000e-9]
    if pick_res[0]:
        index = pick_res[1]
        res_E =  [res_E[index]]
        gn    =  [gn[index]]    
        gg    =  [gg[index]]    
        gfa   =  [gfa[index]]   
        gfb   =  [gfb[index]]
    comp = 'all'
    ap = 0.948 #[barns]
    A = 238
    alpha = ((A-1)/(A+1))**2
    mod_xs_micro = 20.0 #[barns]
    e_bins = np.logspace(0,2,nbins)
    e_freq = np.zeros(nbins-1)
    scat_freq = np.zeros(nbins-1)
    
    captures = 0
    for i in range(neutrons):
        energy = 1000 
        collisions = 0
        while energy > 1:
            e_freq_step = np.histogram(energy,e_bins)[0]
            microXS_elastic = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,energy,'elastic',comp)
            microXS_capture = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,energy,'capture',comp)
            mod_xs_macro     = mod_xs_micro*mf_ratio
            elastic_xs_macro = microXS_elastic
            capture_xs_macro = microXS_capture
            fuel_xs_macro    = elastic_xs_macro+capture_xs_macro
            total_xs_macro   = mod_xs_macro + fuel_xs_macro
            e_freq += e_freq_step/total_xs_macro
            mf_random = np.random.uniform(0,1)
            if mf_random < fuel_xs_macro/total_xs_macro:     #Fuel reaction
                xsi = np.random.uniform(0,1)
                if xsi < capture_xs_macro/fuel_xs_macro:     #Capture in fuel
                    captures += 1; break
                scat_freq += e_freq_step                     #Scatter in fuel
                energy = (energy-alpha*energy)*np.random.uniform(0,1)+alpha*energy    
            else:                                            #Scatter in moderator
                energy = np.random.uniform(0,1)*energy 
    e_freq = np.concatenate([[0],e_freq])/neutrons
    scat_freq = np.concatenate([[0],scat_freq])/neutrons
    abs_ratio = captures/neutrons
    if pick_res[0]:
        return e_bins,e_freq,abs_ratio,scat_freq
    elif by_component:
        xs_tot = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,reaction='elastic') \
               + xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,reaction='capture')
        psi_frac = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,'elastic','psi')/xs_tot
        chi_frac = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,'elastic','chi')/xs_tot
        pot_frac = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,'elastic','pot')/xs_tot
        comp_wise = [e_freq*psi_frac,e_freq*chi_frac,e_freq*pot_frac]
        return e_bins, e_freq, abs_ratio, scat_freq, comp_wise
    else:
        return e_bins,e_freq,abs_ratio


if __name__ == "__main__":

    temp=1000  #temp=1000
    ax = plt.subplot(2,1,1)
    ax.set_xscale("log", nonposx='clip')
    legend_string = []
    neutrons = int(1e4)
    logemin = 0
    logemax = 2 
    nbins = 75


    eBinsVec = []
    eFreqVec = []


    plot_3c = True

    if plot_3c:
        neutrons = int(1e5); logemin = 0; logemax = 2; nbins = 50; temp = 1000;
        mf_ratio = 1000;
        components = ['Total', 'psi', 'chi', 'pot']
        e_bins,e_freq,abs_ratio,scat_freq, comp_wise = zerod_mc(mf_ratio, temp,
                         neutrons, logemin, logemax, nbins, by_component=True)
        ax = plt.subplot(1,1,1)
        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')
        plt.step(e_bins, scat_freq)
        for component in comp_wise:
            plt.step(e_bins, component)
        plt.xlabel('Energy (eV)')
        plt.ylabel('Scatters (#/src neutron)')
        plt.legend(components)
        plt.title('T='+str(temp)+' K  M/F Ratio='+str(mf_ratio))
        plt.show()
    
        # Plot the the components of the xs
        elastic_xs_psi = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,
                                     reaction='elastic',comp='psi')
        elastic_xs_chi = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,
                                     reaction='elastic',comp='chi')
        elastic_xs_pot = xs_from_res(ap,A,res_E,gn,gg,gfa,gfb,temp,e_bins,
                                     reaction='elastic',comp='pot')
        
        ax = plt.subplot(1,1,1)
        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')
        plt.step(e_bins, elastic_xs_pot)
        plt.step(e_bins, elastic_xs_psi)
        plt.step(e_bins, elastic_xs_chi)
        plt.legend(['pot','psi','chi'])
        plt.title('Scattering Components')
        plt.xlabel('eV')
        plt.ylabel('barns')
        plt.show()
        exit()








    for mf_ratio in [10,1e3,1e6]:
        e_bins1,e_freq1,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin,\
                                             logemax, nbins)
        plt.step(e_bins1,e_freq1)
        legend_string.append('mod/fuel ratio: '+str(mf_ratio))
        eBinsVec.append(e_bins1)
        eFreqVec.append(e_freq1)
        print(mf_ratio,"normalized absorption: "+str(abs_ratio))


    plt.legend(legend_string)
    plt.title('T='+str(temp)+' K')

    ax = plt.subplot(2,1,2)
    ax.set_xscale("log")#, nonposx='clip')
    ax.set_yscale("log")#, nonposx='clip')
    plt.step(eBinsVec[0],eFreqVec[0])
    plt.step(eBinsVec[1],eFreqVec[1])
    plt.step(eBinsVec[2],eFreqVec[2])
    plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.xlabel('Energy (eV)')

    plt.show()




