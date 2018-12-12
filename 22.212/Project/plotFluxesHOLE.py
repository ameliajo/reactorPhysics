import matplotlib.pyplot as plt
E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
E_midpoints = [(E_bounds[i]+E_bounds[i+1])/2.0 for i in range(10)]
E_diffs = [(E_bounds[i+1]-E_bounds[i]) for i in range(10)]

def normalize(vec):
	invSum = 1.0/max(vec)
	return [x*invSum for x in vec]

def calcReactionRate(E_diffs,flux,XS):
	RR_denom = sum([E_diffs[i]*flux[i] for i in range(10)])
	RR_numer = sum([XS[i]*flux[i]*E_diffs[i] for i in range(10)])
	return RR_numer/RR_denom



MC_HOLE_k = 1.16307
MC_HOLE_flux = [0.0533311, 0.04768511, 0.02033356, 0.02020343, 0.04833922, 0.02116975, 0.03464084, 0.14510904, 0.30826962, 0.28563944]
MC_HOLE_flux = normalize(MC_HOLE_flux)
MC_HOLE_total = [0.26452482, 0.41823977, 0.54356147, 0.61401882, 0.62318825, 0.42897284, 0.52542178, 0.63591605, 0.92001715, 1.72803998][::-1]
MC_HOLE_absor = [0.01081963, 0.00699548, 0.0528739, 0.17058429, 0.2400336, 0.0537891, 0.14306749, 0.2518359, 0.53323695, 1.33533239][::-1]
MC_HOLE_RR_T = calcReactionRate(E_diffs,MC_HOLE_flux,MC_HOLE_total)
MC_HOLE_RR_A = calcReactionRate(E_diffs,MC_HOLE_flux,MC_HOLE_absor)

MOC_HOLE_k = 1.13031151054
MOC_HOLE_flux = [7.59218404, 8.4304638, 3.87721727, 0.90298816, 0.53355901, 1.23168821, 0.50037014, 0.49255577, 1.07752841, 1.16002763][::-1]
MOC_HOLE_flux = normalize(MOC_HOLE_flux)
MOC_HOLE_total = [0.26452482, 0.41823977, 0.54356147, 0.61401882, 0.62318825, 0.42897284, 0.52542178, 0.63591605, 0.92001715, 1.72803998][::-1]
MOC_HOLE_absor = [0.01081963, 0.00699548, 0.0528739, 0.17058429, 0.2400336, 0.0537891, 0.14306749, 0.2518359, 0.53323695, 1.33533239][::-1]
MOC_HOLE_RR_T = calcReactionRate(E_diffs,MOC_HOLE_flux,MOC_HOLE_total)
MOC_HOLE_RR_A = calcReactionRate(E_diffs,MOC_HOLE_flux,MOC_HOLE_absor)

MOC_TONES_HOLE_k = 1.59946263837
MOC_TONES_HOLE_flux = MOC_fuelFlux0 = [17.0377587, 18.59519274, 8.5210343, 2.00253034, 1.18756419, 2.7313913, 1.11284028, 1.08259866, 2.36731197, 1.48637233][::-gi1]
MOC_TONES_HOLE_flux = normalize(MOC_TONES_HOLE_flux)
TONES_HOLE_total = [0.26408429, 0.41934967, 0.54239907, 0.60479616, 0.61226767, 0.42569651, 0.51841298, 0.60496111, 0.69900817, 2.70680737][::-1]
TONES_HOLE_absor = [0.01007257, 0.00568054, 0.03899287, 0.1290625, 0.20407664, 0.01942985, 0.03307155, 0.05638929, 0.07162838, 0.41157611][::-1]
MOC_TONES_HOLE_RR_T = calcReactionRate(E_diffs,MOC_TONES_HOLE_flux,TONES_HOLE_total)
MOC_TONES_HOLE_RR_A = calcReactionRate(E_diffs,MOC_TONES_HOLE_flux,TONES_HOLE_absor)


print(MC_HOLE_RR_T,MOC_HOLE_RR_T,MOC_HOLE_RR_T)
print(MC_HOLE_RR_A,MOC_HOLE_RR_A,MOC_HOLE_RR_A)





E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
plt.step(E_bounds,[MC_HOLE_flux[0]]+MC_HOLE_flux,label='MC')
plt.step(E_bounds,[MOC_HOLE_flux[0]]+MOC_HOLE_flux,label='MOC using MC-generated XS')
plt.step(E_bounds,[MOC_TONES_HOLE_flux[0]]+MOC_TONES_HOLE_flux,label='MOC using Tones XS')
plt.ylabel('Scalar Flux (normalized)')
plt.xlabel('Energy (eV)')
plt.xscale('log')
plt.title('Neutron flux in 12% enr. fuel pin, with Gd using different methods')
plt.legend(loc='best')
plt.show()














