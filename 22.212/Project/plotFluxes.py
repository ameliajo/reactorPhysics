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

MC_flux = [0.03847132, 0.03439647, 0.01607827, 0.01763081, 0.0445492, 0.01992007, 0.03361357, 0.14907025, 0.32105566, 0.29663482]
MC_flux = normalize(MC_flux)
MC_total = [0.26478544, 0.41821932, 0.56141494, 0.68498781, 0.67487499, 0.48041564, 0.68013205, 0.85794872, 1.07142635, 1.76805619][::-1]
MC_absor = [0.01168467, 0.00893124, 0.07430822, 0.23759845, 0.28911195, 0.10094141, 0.2921717, 0.4682297, 0.67964114, 1.37114701][::-1]
MC_RR_T = calcReactionRate(E_diffs,MC_flux,MC_total)
MC_RR_A = calcReactionRate(E_diffs,MC_flux,MC_absor)


MOC_flux = [13.56930373, 14.40632005, 6.3872181, 1.403151, 0.84202675, 1.89571062, 0.73176022, 0.68137302, 1.48970192, 1.70304391][::-1]
MOC_flux = normalize(MOC_flux)
MOC_total = [0.26478544, 0.41821932, 0.56141494, 0.68498781, 0.67487499, 0.48041564, 0.68013205, 0.85794872, 1.07142635, 1.76805619][::-1]
MOC_absor = [0.01168467, 0.00893124, 0.07430822, 0.23759845, 0.28911195, 0.10094141, 0.2921717, 0.4682297, 0.67964114, 1.37114701][::-1]
MOC_RR_T = calcReactionRate(E_diffs,MOC_flux,MOC_total)
MOC_RR_A = calcReactionRate(E_diffs,MOC_flux,MOC_absor)

MOC_tones_flux = [27.61267839, 29.66454511, 13.42243859, 3.02199575, 1.82631201, 4.06239555, 1.57222825, 1.39002337, 2.61078589, 1.12811979][::-1]
MOC_tones_flux = normalize(MOC_tones_flux)
tones_total  = [0.26402023, 0.42107671, 0.55836669, 0.65528308, 0.64502773, 0.47569187, 0.66881908, 0.86016566, 1.05419221, 5.04548983][::-1]
tones_absor  = [0.00959068, 0.00600729, 0.0435695, 0.14318842, 0.21172985, 0.02872196, 0.05342193, 0.09899989, 0.11983258, 0.77415291][::-1]
MOC_tones_RR_T = calcReactionRate(E_diffs,MOC_tones_flux,tones_total)
MOC_tones_RR_A = calcReactionRate(E_diffs,MOC_tones_flux,tones_absor)



E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
plt.step(E_bounds,[MC_flux[0]]+MC_flux,label='MC')
plt.step(E_bounds,[MOC_flux[0]]+MOC_flux,label='MOC using MC-generated XS')
plt.step(E_bounds,[MOC_tones_flux[0]]+MOC_tones_flux,label='MOC using Tones XS')
plt.ylabel('Scalar Flux (normalized)')
plt.xlabel('Energy (eV)')
plt.xscale('log')
plt.title('Neutron flux in 9% enr. fuel pin using different methods')
plt.legend(loc='best')
plt.show()













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


MC_k = 1.60629
MC_flux = [0.02879449, 0.02778762, 0.0140252, 0.01609068, 0.04363773, 0.01945098, 0.03246101, 0.14740258, 0.32073465, 0.29916951]
MC_flux = normalize(MC_flux)
MC_total = [0.26482608, 0.41782873, 0.5797005, 0.72601689, 0.70774452, 0.51290634, 0.77004306, 1.00974167, 1.28751244, 2.18380972][::-1]
MC_absor = [0.01209566, 0.01013866, 0.08950889, 0.27487643, 0.32051569, 0.13075438, 0.37871566, 0.61644707, 0.89190662, 1.78308371][::-1]
MC_RR_T = calcReactionRate(E_diffs,MC_flux,MC_total)
MC_RR_A = calcReactionRate(E_diffs,MC_flux,MC_absor)

MOC_k = 1.60706049226
MOC_flux = [16.94125253, 17.8739986, 7.80542238, 1.68968712, 1.00322899, 2.26142942, 0.82319294, 0.72828042, 1.48691383, 1.57564102][::-1]
MOC_flux = normalize(MOC_flux)
MOC_total = [0.26482608, 0.41782873, 0.5797005, 0.72601689, 0.70774452, 0.51290634, 0.77004306, 1.00974167, 1.28751244, 2.18380972][::-1]
MOC_absor = [0.01209566, 0.01013866, 0.08950889, 0.27487643, 0.32051569, 0.13075438, 0.37871566, 0.61644707, 0.89190662, 1.78308371][::-1]
MOC_RR_T = calcReactionRate(E_diffs,MOC_flux,MOC_total)
MOC_RR_A = calcReactionRate(E_diffs,MOC_flux,MOC_absor)

MOC_tones_k = 1.7629326
MOC_tones_flux = [34.04768914, 36.68243303, 16.61561472, 3.73428848, 2.2341029, 4.97787129, 1.82204522, 1.53810266, 2.77921125, 1.03687844][::-1]
MOC_tones_flux = normalize(MOC_tones_flux)
tones_total  = [0.26379803, 0.42145659, 0.57197067, 0.68110208, 0.66157794, 0.5003477, 0.7547992, 1.01321759, 1.26616191, 6.41545478][::-1]
tones_absor  = [0.00928902, 0.00620524, 0.04754651, 0.15140082, 0.2179899, 0.0332696, 0.06470017, 0.12455949, 0.14860541, 0.98589298][::-1]
MOC_tones_RR_T = calcReactionRate(E_diffs,MOC_tones_flux,tones_total)
MOC_tones_RR_A = calcReactionRate(E_diffs,MOC_tones_flux,tones_absor)


# print(MC_RR_T,MOC_RR_T,MOC_tones_RR_T)
# print(MC_RR_A,MOC_RR_A,MOC_tones_RR_A)

# E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
# plt.step(E_bounds,[MC_flux[0]]+MC_flux,label='MC')
# plt.step(E_bounds,[MOC_flux[0]]+MOC_flux,label='MOC using MC-generated XS')
# plt.step(E_bounds,[MOC_tones_flux[0]]+MOC_tones_flux,label='MOC using Tones XS')
# # plt.step(E_bounds,[tonesTheory[0]]+tonesTheory,label='Tones Eq. Solution')
# plt.ylabel('Scalar Flux (normalized)')
# plt.xlabel('Energy (eV)')
# plt.xscale('log')
# plt.title('Neutron flux in 12% enr. fuel pin using different methods')
# plt.legend(loc='best')
# plt.show()









