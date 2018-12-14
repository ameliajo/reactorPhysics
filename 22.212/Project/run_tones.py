from math import pi
import matplotlib.pyplot as plt
from numpy import ma
import sys
sys.path.append('./dilutionTables/')
sys.path.append('./GettingXS/')
sys.path.append('./collisionProb/')
from interpret import *
from import_XS_help import *
from XS_nuclideSpecific import *
from getCollisionProb import *
from XS import modTotal0
import matplotlib.colors as colors
import matplotlib.cm as cmx

class Pin:
    def __init__(self,N,radius,all_nuclides_pin):
        self.U235 = Nuclide(N['U235'],radius['U235'],all_nuclides_pin[0])
        self.U238 = Nuclide(N['U238'],radius['U238'],all_nuclides_pin[1])
        self.O16  = Nuclide(N['O16' ],radius['O16' ],all_nuclides_pin[2])
        self.nuclides = [self.U235,self.U238,self.O16]

    def calcSigT(self,nGroups):
        self.SigT = []
        for g in range(nGroups):
            self.SigT.append(self.U235.SigT[g]        + \
                             self.U238.openMC_SigT[g] + \
                             self.O16.openMC_SigT[g])




def run_tones():

        #jet = cm = plt.get_cmap('hot')
        jet = cm = plt.get_cmap('tab10')
        #jet = cm = plt.get_cmap('autumn')
        cNorm  = colors.Normalize(vmin=0, vmax=10)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

        # Define problem geometry
        pinRad = 0.39128;  pitch = 1.26
        l_bar  = 2*pinRad; C     = 0.14


        E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
        E_bounds = [1e-5,4,10,18,25,34,45,5530,821e3,2.e7]
        E_bounds_inc_g = E_bounds[::-1]

        resonanceGroups = [1,2,3,4,5,6,7]
        nonResonanceGroups = [0,8]


        ###############################################################################
        # Load in dilution table
        ###############################################################################

        # groupsU235[g].sigT[d] will give you sigT for group g for dilution index d
        grouprDataU8 = interpretGENDF('u238',False,'dilutionTables/u238_tape26')
        grouprDataU8.reverse() # Because NJOY does groups in inverse order (low -> high)
        dilution = grouprDataU8[0].dilutionVals # All dilution values are the same, it
                                              # doesn't matter I'm pulling this from 
                                              # energy group 0
        nGroups = len(grouprDataU8)

        # Gotten from OpenMC input, using uo2.get_nuclide_atom_densities and multiplying
        # by 1E24 bc the number densities are provided in atoms / b cm
        N_dict_Hi = { 'U235' : N_hi_U235, 'U238' : N_hi_U238, 'O16'  : N_hi_O16 }
        N_dict_Lo = { 'U235' : N_lo_U235, 'U238' : N_lo_U238, 'O16'  : N_lo_O16 }

        # Called ``AP'' in ENDF, found in MF 2 MT 151
        radius = { 'U234': 8.930000E-1, 'U235': 9.602000E-1, \
                   'U238': 9.480000E-1, 'U236': 9.354000E-1, \
                   'O16' : 5.562563E-1, 'O17' : 5.780000E-1 }


        # Define pins materials
        # Load in openMC data in for all materials. 
        pins = [Pin([N_dict_Hi,N_dict_Lo][i%2],radius,all_nuclides[i]) for i in range(9)]

        
        """
        for g in resonanceGroups:
            plt.plot(grouprDataU8[g].dilutionVals,grouprDataU8[g].sigT,label=str(int(grouprDataU8[g].E_low))+" - "+str(int(grouprDataU8[g].E_high))+' eV',color=scalarMap.to_rgba(g))

        plt.xscale('log')
        plt.xlabel('background XS (b)')
        plt.ylabel('absorption XS (b)')
        plt.ylabel('total XS (b)')
        plt.legend(loc='best')
        plt.title('U-238 sigT dependency on background XS')
        plt.show()
        return
        """



        ###############################################################################
        # Assume initial background cross sections for resonance nuclides.  
        # They can be evaluated using the conventional equivalence method with the 
        # Dancoff  correction.
        ###############################################################################

        #------------------------------------------------------------------------------
        # We're going to treat U-238 as the resonant nuclide for now, for pin0
        #------------------------------------------------------------------------------
        nonRes = [pins[0].U235,pins[0].O16]

        # background = sum_nonRes [N_nonRes * sig_potNonRes / N_res] + 1/(N_res*l_bar)
        #sig0 = sum([nuclide.N*nuclide.pot for nuclide in nonRes])/pins[0].U238.N + \
        #       1.0 / (pins[0].U238.N*l_bar*(1.0-C))
        sig0 = sum([nuclide.N*nuclide.pot for nuclide in nonRes])/pins[0].U238.N + \
               1.0 / (pins[0].U238.N*l_bar)



        sig0EnergyVec = [sig0]*nGroups

        converged = False
        counter = 0

            
        absorptionXS_g0 = []
        absorptionXS_g1 = []
        absorptionXS_g2 = []
        absorptionXS_g3 = []
        absorptionXS_g4 = []
        absorptionXS_g5 = []
        absorptionXS_g6 = []
        absorptionXS_g7 = []
        absorptionXS_g8 = []


        while not converged:

            ###########################################################################
            # Evaluate the effective cross sections of resonance nuclides using the 
            # conventional equivalence theory.
            ###########################################################################

            boundingDilutionsVec = [None]*nGroups
            #for g in range(nGroups):
            # Only look up the dilution data for the groups that are in the 
            # resonance region
            for g in resonanceGroups:
                for i in range(len(dilution)-1):
                    if dilution[i] > sig0EnergyVec[g] > dilution[i+1]:
                        boundingDilutionsVec[g] = [(i,dilution[i]),(i+1,dilution[i+1])]
                        break


            if None in [boundingDilutionsVec[g] for g in resonanceGroups]: 
                raise ValueError('Ideal dilution value out of range')

            # Pulling cross sections from the dilution table
            for pin in pins:
                pin.U238.resetXS(nGroups)
                for g in resonanceGroups:
                    # Returns sigT, sigF, sigA, nuBar
                    pin.U238.addXS(getDataFromEqTable(boundingDilutionsVec[g],\
                                   grouprDataU8[g],sig0EnergyVec[g]),g)

                # Making these microscopic cross sections into macroscropic cross sections
                pin.U238.convertToMacro(nGroups)

            
            absorptionXS_g0.append(pins[0].U238.SigA[0])
            absorptionXS_g1.append(pins[0].U238.SigA[1])
            absorptionXS_g2.append(pins[0].U238.SigA[2])
            absorptionXS_g3.append(pins[0].U238.SigA[3])
            absorptionXS_g4.append(pins[0].U238.SigA[4])
            absorptionXS_g5.append(pins[0].U238.SigA[5])
            absorptionXS_g6.append(pins[0].U238.SigA[6])
            absorptionXS_g7.append(pins[0].U238.SigA[7])
            absorptionXS_g8.append(pins[0].U238.SigA[8])

            

            ###########################################################################
            # Evaluate group-wise collision probability using the effctive XS from dilution
            ###########################################################################

            #--------------------------------------------------------------------------
            # Calculate the macro. SigT for the high/low enr. fuel pins, and plug into the
            # collision probability MC script
            #--------------------------------------------------------------------------
            SigT_hi = [ sum(nucl.SigT[g] for nucl in pins[0].nuclides) for g in range(nGroups) ]
            SigT_lo = [ sum(nucl.SigT[g] for nucl in pins[1].nuclides) for g in range(nGroups) ]

            # For hi/lo enr. of fuel, use values we just pulled from dilution table
            # For moderator, use openMC values generated from grid_3x3.py
            collisionProbs = [                                                   \
                getCollisionProb( pitch, pinRad, plot=False, numParticles=5000,  \
                  hole=True, fSigT_hi=SigT_hi[g], fSigT_lo=SigT_lo[g],          \
                  mSigT=modTotal0[g], verbose=False, startNeutronsFrom=0 )       \
                for g in range(nGroups)]





            ###########################################################################
            # Update the background cross section 
            ###########################################################################

            # sig0 = SUM_pins SUM_nonRes P(pin->me) * V(pin) * N(nonResInPin)*sig(pot)
            #        ------------------------------------------------------------------
            #                 SUM_pins * P(pin->me) * V(pin) * N(resInPin)

            tones_Numer = 0.0
            tones_Denom = 0.0

            for pinID,pin in enumerate(pins):
                pin.calcSigT(nGroups)
                P_0_to_i = np.array([collisionProbs[g][pinID] for g in range(nGroups)])
                # RECIPROCITY RELATION
                # P(i->0) = P(0->i) * SigT_0 / SigT_i
                P_i_to_0 = P_0_to_i * pins[0].SigT / pin.SigT
                SUM_nonRes_sigPot = pin.U235.N * pin.U235.pot * 1E-24 + \
                                    pin.O16.N  * pin.O16.pot  * 1E-24
                
                tones_Numer += P_i_to_0 * SUM_nonRes_sigPot
                tones_Denom += P_i_to_0 * pin.U238.N * 1e-24


            #plt.step(E_bounds,[1e24*tones_Numer[0]/tones_Denom[0]]+list(1e24*tones_Numer/tones_Denom),label='inter'+str(counter),color=scalarMap.to_rgba(counter))

            sig0EnergyVec = tones_Numer/tones_Denom
            ###########################################################################
            # Repeat
            ###########################################################################

            counter += 1
            if counter > 3:
                print(counter)
                plt.plot(absorptionXS_g0)
                plt.plot(absorptionXS_g1)
                plt.plot(absorptionXS_g2)
                plt.plot(absorptionXS_g3)
                plt.plot(absorptionXS_g4)
                plt.plot(absorptionXS_g5)
                plt.plot(absorptionXS_g6)
                plt.plot(absorptionXS_g7)
                plt.plot(absorptionXS_g8)
                SigT_hi = [ sum(nucl.SigT[g] for nucl in pins[0].nuclides) for g in range(nGroups) ]
                SigT_lo = [ sum(nucl.SigT[g] for nucl in pins[1].nuclides) for g in range(nGroups) ]
                SigA_hi = [ sum(nucl.SigA[g] for nucl in pins[0].nuclides) for g in range(nGroups) ]
                SigA_lo = [ sum(nucl.SigA[g] for nucl in pins[1].nuclides) for g in range(nGroups) ]
                nuSigF_hi = [ sum(nucl.nuSigF[g] for nucl in pins[0].nuclides) for g in range(nGroups) ]
                nuSigF_lo = [ sum(nucl.nuSigF[g] for nucl in pins[1].nuclides) for g in range(nGroups) ]
                plt.show()

                return


                """
                f = open("XS-tones.py","a")
                for i in range(9):
                    f.write("fuelTotal"+str(i)+"  = "+str([float("%.8f"%[SigT_hi,SigT_lo][i%2][g]) for g in range(nGroups)])+"\n")
                    f.write("fuelAbsorption"+str(i)+"  = "+str([float("%.8f"%[SigA_hi,SigA_lo][i%2][g]) for g in range(nGroups)])+"\n")
                    f.write("fuelNuFission"+str(i)+"  = "+str([float("%.8f"%[nuSigF_hi,nuSigF_lo][i%2][g]) for g in range(nGroups)])+"\n")

                    f.write("\n\n")
                f.close()

                frac = []
                for i in range(10):
                    numer = sig0EnergyVec[i]+pins[0].U235.pot
                    denom = sig0EnergyVec[i]+pins[0].U235.sigT[i]
                    frac.append(numer/denom)

                print(frac)

                break
                """

        """
        ax = plt.gca()
        ax.set_facecolor('xkcd:pale grey')
        ax.set_facecolor('xkcd:off white')
        plt.xscale('log')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Sig0 approximation')
        plt.legend(loc='best')
        #plt.savefig('sig0Estimations.png')
        plt.show()

        """




if __name__ == "__main__":
    run_tones()

