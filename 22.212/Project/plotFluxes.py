import matplotlib.pyplot as plt
"""

########################################################################
# 6 % Enrichment in high pin. 4 % in low pin
########################################################################
# MC
# mc = [0.05627205, 0.04522832, 0.01901171, 0.01913341, 0.04762022, 0.02118881, 0.03418096, 0.15094522, 0.318804, 0.29161574]
MC_fuelFlux0 = [0.055448, 0.04487209, 0.01875088, 0.01928183, 0.04722538, 0.02137325, 0.03567234, 0.15209708, 0.32118123, 0.29461401]

sum_MC_fuelFlux0 = sum(MC_fuelFlux0)
MC_fuelFlux0 = [x / sum_MC_fuelFlux0 for x in MC_fuelFlux0]
# MC_fuelFlux0.reverse()

MOC_fuelFlux0 = [10.09771892, 10.89388971, 4.93458145, 1.13330597, 0.66698134, 1.53610706, 0.6059869, 0.60868883, 1.50008976, 1.86237245]
sum_MOC_fuelFlux0 = sum(MOC_fuelFlux0)
MOC_fuelFlux0 = [x / sum_MOC_fuelFlux0 for x in MOC_fuelFlux0]
MOC_fuelFlux0.reverse()


tonesTheory = [ 1.01394269, 0.99758961, 0.93844619, 0.84925338, 0.8755602, 0.87954041, 0.70281634, 0.58765548, 0.50331242, 0.11789354]
sum_tonesTheory = sum(tonesTheory)
tonesTheory.reverse()
tonesTheory = [x / sum_tonesTheory for x in tonesTheory]
E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
for i in range(10):
	tonesTheory[i] /= (E_bounds[i+1]+E_bounds[i])
tonesTheory.reverse()



tonesMOC = [21.05819301, 22.77673383, 10.39519434, 2.42790356, 1.43375061, 3.27890329, 1.29671076, 1.21592632, 2.46986688, 1.2869658 ]
sum_tonesMOC = sum(tonesMOC)
tonesMOC = [x / sum_tonesMOC for x in tonesMOC]

print(MC_fuelFlux0)
print(MOC_fuelFlux0)
print(tonesMOC)
print(tonesTheory)
tonesMOC.reverse()

# MOC_fuelFlux0.reverse()
# plt.plot(MC_fuelFlux0,label='MC')
# plt.plot(MOC_fuelFlux0,label='MOC using MC-generated XS')
# plt.plot(tonesMOC,label='MOC using Tones XS')
# plt.plot(tonesTheory,label='Tones Eq. Solution')
# plt.title('Neutron flux in 6% enr. fuel pin using different methods')
# plt.xlabel('Energy (eV)')
# plt.ylabel('Scalar Flux (normalized)')
# plt.legend(loc='best')
# plt.show()

E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
plt.step(E_bounds,[MC_fuelFlux0[0]]+MC_fuelFlux0,label='MC')
plt.step(E_bounds,[MOC_fuelFlux0[0]]+MOC_fuelFlux0,label='MOC using MC-generated XS')
plt.step(E_bounds,[tonesMOC[0]]+tonesMOC,label='MOC using Tones XS')
plt.step(E_bounds,[tonesTheory[0]]+tonesTheory,label='Tones Eq. Solution')
plt.ylabel('Scalar Flux (normalized)')
plt.xlabel('Energy (eV)')
plt.xscale('log')
plt.title('Neutron flux in 6% enr. fuel pin using different methods')
plt.legend(loc='best')
plt.show()


"""


########################################################################
# 9 % Enrichment in high pin. 4 % in low pin
########################################################################
# MC
MC_fuelFlux0 = [0.03824895, 0.03413355, 0.01628323, 0.01757013, 0.04479519, 0.01995093, 0.03380981, 0.14810622, 0.32018227, 0.29649381]
sum_MC_fuelFlux0 = sum(MC_fuelFlux0)
MC_fuelFlux0 = [x / sum_MC_fuelFlux0 for x in MC_fuelFlux0]
# MC_fuelFlux0.reverse()

MOC_fuelFlux0 = [13.5196022, 14.40633473, 6.46404503, 1.42415558, 0.84707415, 1.90625945, 0.73329303, 0.68103082, 1.51879883, 1.71886332]
sum_MOC_fuelFlux0 = sum(MOC_fuelFlux0)
MOC_fuelFlux0 = [x / sum_MOC_fuelFlux0 for x in MOC_fuelFlux0]
MOC_fuelFlux0.reverse()


tonesTheory = [ 1.01999297, 0.99650289, 0.91326466, 0.80437664, 0.83822861, 0.83392992, 0.61710356, 0.49083038, 0.40640608, 0.08524643]
sum_tonesTheory = sum(tonesTheory)
tonesTheory = [x / sum_tonesTheory for x in tonesTheory]
tonesTheory.reverse()
E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
for i in range(10):
	tonesTheory[i] /= (E_bounds[i+1]+E_bounds[i])
tonesTheory.reverse()


tonesMOC = [27.54297667, 29.67489576, 13.500035, 3.04042113, 1.8155832, 4.03846653, 1.55388277, 1.37043008, 2.63352521, 1.1251432]
sum_tonesMOC = sum(tonesMOC)
tonesMOC = [x / sum_tonesMOC for x in tonesMOC]
tonesMOC.reverse()

# MOC_fuelFlux0.reverse()
# plt.plot(MC_fuelFlux0,label='MC')
# plt.plot(MOC_fuelFlux0,label='MOC using MC-generated XS')
# plt.plot(tonesMOC,label='MOC using Tones XS')
# plt.plot(tonesTheory,label='Tones Eq. Solution')
# plt.title('Neutron flux in 6% enr. fuel pin using different methods')
# plt.xlabel('Energy (eV)')
# plt.ylabel('Scalar Flux (normalized)')
# plt.legend(loc='best')
# plt.show()

E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
plt.step(E_bounds,[MC_fuelFlux0[0]]+MC_fuelFlux0,label='MC')
plt.step(E_bounds,[MOC_fuelFlux0[0]]+MOC_fuelFlux0,label='MOC using MC-generated XS')
plt.step(E_bounds,[tonesMOC[0]]+tonesMOC,label='MOC using Tones XS')
plt.step(E_bounds,[tonesTheory[0]]+tonesTheory,label='Tones Eq. Solution')
plt.ylabel('Scalar Flux (normalized)')
plt.xlabel('Energy (eV)')
plt.xscale('log')
plt.title('Neutron flux in 9% enr. fuel pin using different methods')
plt.legend(loc='best')
plt.show()











