import numpy as np
import matplotlib.pyplot as plt
from math import log

# My file imports here
from res_construction import xs_from_res
from nslowing_resonances import zerod_mc
from slbw import full_xs
from nslowing import nslowing

#Resonance data
res_E = [6.673491e+0, 2.087152e+1, 3.668212e+1]
J = [5.000000e-1, 5.000000e-1, 5.000000e-1]
gn = [1.475792e-3, 1.009376e-2, 3.354568e-2]
gg = [2.300000e-2, 2.286379e-2, 2.300225e-2]
gfa = [0.000000, 5.420000e-8, 0.000000]
gfb = [9.990000e-9, 0.000000, 9.770000e-9]
ap = 0.948 #[barns]
A = 238


plot_1a = False
plot_2abc = False
plot_3a = False
plot_3b = False
plot_3c = True

if plot_1a:
    points = 100000
    energies = np.linspace(0.01,100,points); temp = 1000;
    capture_xs, elastic_xs, total_xs = full_xs(energies,temp)
    plt.loglog(energies,capture_xs,'C0', linewidth=5.0)
    plt.xlabel('Energy [eV]')
    plt.ylabel('Cross Section [barns]')
    plt.title(str(temp)+' K Cross Sections')
    plt.show()
    
    
    #Plot single cross section to check
    # points = 100000
    # energies = np.logspace(-1,3,points)
    # capture_xs = np.zeros([points])
    # elastic_xs = np.zeros([points])
    # total_xs = np.zeros([points])
    # temp = 0
    # for j in range(points):
    #     capture_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
    #                                      gfa,gfb,temp,energies[j],
    #                                      reaction='capture')
    #     elastic_xs[j]=xs_from_res(ap,A,res_E,J,gn,gg,
    #                                      gfa,gfb,temp,energies[j],
    #                                      reaction='elastic')
    # total_xs =capture_xs+elastic_xs
    # plt.loglog(energies,elastic_xs)
    # plt.show()

if plot_2abc:
    bins = 100
    e_bins, leth_bins, e_freq, leth_freq, collision_dist = nslowing(bins)
    ax1 = plt.subplot(3,1,1)
    plt.step(e_bins,e_freq)
    plt.title('Equal Energy Bins')
    ax2 = plt.subplot(3,1,2)
    ax2.set_xscale("log", nonposx='clip')
    plt.step(leth_bins,leth_freq)
    plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.title('Equal Lethargy Bins')
    ax3=plt.subplot(3,1,3)
    ax3.set_xscale("log", nonposx='clip')
    ax3.set_yscale("log", nonposy='clip')
    for i in range(3):
        plt.step(leth_bins,collision_dist[i+1])
    plt.legend(['1 Collison','2 Collisions','3 Collisions'])
    plt.title('Energy Distribution after collisions (Equal Lethargy)')
    plt.xlabel('Energy [eV]')
    plt.show()

if plot_3a or plot_3b:
    if plot_3a:
        temp=0
    else:
        temp=1000
    ax = plt.subplot(2,1,1)
    ax.set_xscale("log", nonposx='clip')
    legend_string = []
    neutrons = int(1e4); logemin = 0; logemax = 2; nbins = 75;
    mf_ratio = 10; 
    e_bins1,e_freq1,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin,
                                         logemax, nbins)
    plt.step(e_bins1,e_freq1)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+ \
                         ', Normalized Absorption='+str(abs_ratio))
    mf_ratio = 1000;
    e_bins2,e_freq2,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin,
                                         logemax, nbins)
    plt.step(e_bins2,e_freq2)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+ \
                         ', Normalized Absorption='+str(abs_ratio))
    mf_ratio = 1000000;
    e_bins3,e_freq3,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin,
                                         logemax, nbins)
    plt.step(e_bins3,e_freq3)
    legend_string.append('Mod/Fuel Ratio: '+str(mf_ratio)+ \
                         ', Normalized Absorption='+str(abs_ratio))
    plt.legend(legend_string)
    plt.title('T='+str(temp)+' K')

    ax = plt.subplot(2,1,2)
    ax.set_xscale("log")#, nonposx='clip')
    ax.set_yscale("log")#, nonposx='clip')
    plt.step(e_bins1,e_freq1)
    plt.step(e_bins2,e_freq2)
    plt.step(e_bins3,e_freq3)
    plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.xlabel('Energy (eV)')

    plt.show()

if plot_3c:
    neutrons = 100000; logemin = 0; logemax = 2; nbins = 75; temp = 1000;
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

# e_bin_width = np.zeros([len(e_bins)])
# for i in range(len(e_bins)-2):
#     e_bin_width[i+1] = e_bins[i+2]-e_bins[i+1]
# 
# ax.set_xscale("log", nonposx='clip')
# plt.step(e_bins,e_freq/e_bin_width)
# plt.ylabel('Normalized Flux/eV')
# plt.show()

