import numpy as np
import matplotlib.pyplot as plt
from math import log
#from res_construction import xs_from_res
from nslowing_resonances import zerod_mc
#from slbw import full_xs
#from nslowing import nslowing




if __name__ == "__main__":

    temp=0
    temp=1000
    ax = plt.subplot(2,1,1)
    ax.set_xscale("log", nonposx='clip')
    legend_string = []
    neutrons = int(1e4)
    logemin = 0
    logemax = 2 
    nbins = 75


    eBinsVec = []
    eFreqVec = []


    for mf_ratio in [10,1e3,1e6]:
        e_bins1,e_freq1,abs_ratio = zerod_mc(mf_ratio, temp, neutrons, logemin,\
                                             logemax, nbins)
        plt.step(e_bins1,e_freq1)
        legend_string.append('mod/fuel ratio: '+str(mf_ratio))#+ \
        #                     ', normalized absorption='+str(abs_ratio))
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




