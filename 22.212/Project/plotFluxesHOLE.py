
# The only difference between the odd and even pins were that the even ones had gd. both had enrichment of 4%

# import matplotlib.pyplot as plt
# E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
# E_midpoints = [(E_bounds[i]+E_bounds[i+1])/2.0 for i in range(10)]
# E_diffs = [(E_bounds[i+1]-E_bounds[i]) for i in range(10)]

# def normalize(vec):
# 	invSum = 1.0/max(vec)
# 	return [x*invSum for x in vec]

# def calcReactionRate(E_diffs,flux,XS):
# 	RR_denom = sum([E_diffs[i]*flux[i] for i in range(10)])
# 	RR_numer = sum([XS[i]*flux[i]*E_diffs[i] for i in range(10)])
# 	return RR_numer/RR_denom



# MC_HOLE_k = 1.16307
# MC_HOLE_flux = [0.0533311, 0.04768511, 0.02033356, 0.02020343, 0.04833922, 0.02116975, 0.03464084, 0.14510904, 0.30826962, 0.28563944]
# MC_HOLE_flux = normalize(MC_HOLE_flux)
# MC_HOLE_total = [0.26452482, 0.41823977, 0.54356147, 0.61401882, 0.62318825, 0.42897284, 0.52542178, 0.63591605, 0.92001715, 1.72803998][::-1]
# MC_HOLE_absor = [0.01081963, 0.00699548, 0.0528739, 0.17058429, 0.2400336, 0.0537891, 0.14306749, 0.2518359, 0.53323695, 1.33533239][::-1]
# MC_HOLE_RR_T = calcReactionRate(E_diffs,MC_HOLE_flux,MC_HOLE_total)
# MC_HOLE_RR_A = calcReactionRate(E_diffs,MC_HOLE_flux,MC_HOLE_absor)

# MOC_HOLE_k = 1.13031151054
# MOC_HOLE_flux = [7.59218404, 8.4304638, 3.87721727, 0.90298816, 0.53355901, 1.23168821, 0.50037014, 0.49255577, 1.07752841, 1.16002763][::-1]
# MOC_HOLE_flux = normalize(MOC_HOLE_flux)
# MOC_HOLE_total = [0.26452482, 0.41823977, 0.54356147, 0.61401882, 0.62318825, 0.42897284, 0.52542178, 0.63591605, 0.92001715, 1.72803998][::-1]
# MOC_HOLE_absor = [0.01081963, 0.00699548, 0.0528739, 0.17058429, 0.2400336, 0.0537891, 0.14306749, 0.2518359, 0.53323695, 1.33533239][::-1]
# MOC_HOLE_RR_T = calcReactionRate(E_diffs,MOC_HOLE_flux,MOC_HOLE_total)
# MOC_HOLE_RR_A = calcReactionRate(E_diffs,MOC_HOLE_flux,MOC_HOLE_absor)

# MOC_TONES_HOLE_k = 1.59946263837
# MOC_TONES_HOLE_flux = [17.0377587, 18.59519274, 8.5210343, 2.00253034, 1.18756419, 2.7313913, 1.11284028, 1.08259866, 2.36731197, 1.48637233][::-1]
# MOC_TONES_HOLE_flux = normalize(MOC_TONES_HOLE_flux)
# TONES_HOLE_total = [0.26408429, 0.41934967, 0.54239907, 0.60479616, 0.61226767, 0.42569651, 0.51841298, 0.60496111, 0.69900817, 2.70680737][::-1]
# TONES_HOLE_absor = [0.01007257, 0.00568054, 0.03899287, 0.1290625, 0.20407664, 0.01942985, 0.03307155, 0.05638929, 0.07162838, 0.41157611][::-1]
# MOC_TONES_HOLE_RR_T = calcReactionRate(E_diffs,MOC_TONES_HOLE_flux,TONES_HOLE_total)
# MOC_TONES_HOLE_RR_A = calcReactionRate(E_diffs,MOC_TONES_HOLE_flux,TONES_HOLE_absor)


# print(MC_HOLE_RR_T,MOC_HOLE_RR_T,MOC_HOLE_RR_T)
# print(MC_HOLE_RR_A,MOC_HOLE_RR_A,MOC_HOLE_RR_A)





# E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
# plt.step(E_bounds,[MC_HOLE_flux[0]]+MC_HOLE_flux,label='MC')
# plt.step(E_bounds,[MOC_HOLE_flux[0]]+MOC_HOLE_flux,label='MOC using MC-generated XS')
# plt.step(E_bounds,[MOC_TONES_HOLE_flux[0]]+MOC_TONES_HOLE_flux,label='MOC using Tones XS')
# plt.ylabel('Scalar Flux (normalized)')
# plt.xlabel('Energy (eV)')
# plt.xscale('log')
# plt.title('Neutron flux in 4% enr. fuel pin, with Gd using different methods')
# plt.legend(loc='best')
# plt.show()








# 9% enrichment on even pins, 4% on odd pins. Even pins have 0.0007 Gd in them


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



MC_HOLE_k = 1.35213
MC_HOLE_flux = [0.03395539, 0.03445963, 0.01687826, 0.0177109, 0.04460918, 0.02006118, 0.03341757, 0.14425869, 0.31106229, 0.30028612]
MC_HOLE_flux = normalize(MC_HOLE_flux)
MC_HOLE_total = [0.26487685, 0.41772471, 0.56531272, 0.68413602, 0.68806162, 0.48208101, 0.68292823, 0.88993153, 1.27721041, 2.42335511][::-1]
MC_HOLE_absor = [0.01168578, 0.00892975, 0.07616205, 0.23803417, 0.30182539, 0.10258097, 0.29494267, 0.49988722, 0.8840958, 2.02428324][::-1]
MC_HOLE_RR_T = calcReactionRate(E_diffs,MC_HOLE_flux,MC_HOLE_total)
MC_HOLE_RR_A = calcReactionRate(E_diffs,MC_HOLE_flux,MC_HOLE_absor)

MOC_HOLE_k = 1.33934387795
MOC_HOLE_flux = [11.47086666, 13.25787703, 6.59323292, 1.59689152, 0.96773537, 2.00669619, 0.83586748, 0.82828083, 1.88500417, 2.66910848][::-1]
MOC_HOLE_flux = normalize(MOC_HOLE_flux)
MOC_HOLE_total = [0.26487685, 0.41772471, 0.56531272, 0.68413602, 0.68806162, 0.48208101, 0.68292823, 0.88993153, 1.27721041, 2.42335511][::-1]
MOC_HOLE_absor = [0.01168578, 0.00892975, 0.07616205, 0.23803417, 0.30182539, 0.10258097, 0.29494267, 0.49988722, 0.8840958, 2.02428324][::-1]
MOC_HOLE_RR_T = calcReactionRate(E_diffs,MOC_HOLE_flux,MOC_HOLE_total)
MOC_HOLE_RR_A = calcReactionRate(E_diffs,MOC_HOLE_flux,MOC_HOLE_absor)

MOC_TONES_HOLE_k = 1.70633319684
MOC_TONES_HOLE_flux = [27.68604478, 29.86572646, 13.44103488, 3.06133303, 1.80656744, 4.11674109, 1.56887981, 1.38712896, 2.66356837, 1.1293358][::-1]
MOC_TONES_HOLE_flux = normalize(MOC_TONES_HOLE_flux)
TONES_HOLE_total  = [0.26401177, 0.4204555, 0.56116607, 0.65559848, 0.65978078, 0.47555222, 0.66849357, 0.85983972, 1.05362548, 5.04471841][::-1]
TONES_HOLE_absor  = [0.00959169, 0.00599961, 0.04466119, 0.14482695, 0.2260347, 0.02872312, 0.05337423, 0.09898159, 0.11969691, 0.77336633][::-1]
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
plt.title('Neutron flux in 9% enr. fuel pin, with Gd using different methods')
plt.legend(loc='best')
plt.show()


























