
import matplotlib.pyplot as plt

# MC_flux 		= [0.2915480643163238, 0.3178388084844844, 0.15051425851121283, 0.035301110346496314, 0.021150826010103416, 0.04673392187154586, 0.019081170691700716, 0.018555746103953677, 0.04440512174328665, 0.05487097192089243]
# MOC_flux 		= [0.29839839711521937, 0.3219260957418274, 0.1458221611019546, 0.03349040387916729, 0.019710012165971555, 0.04539360702480042, 0.017907561209162745, 0.017987406131318443, 0.04432925727674683, 0.055035098353831385]
# MOC_tones_flux 	= [0.3113268304124537, 0.3367339423223383, 0.15368378967069207, 0.03589441503945606, 0.021196739568359668, 0.048475696277449326, 0.019170726124545286, 0.017976399353967117, 0.03651480575403351, 0.019026655476705016]
# tones_flux 		= [0.27220654333409494, 0.34047393864576864, 0.18740661552288226, 0.1040170872848463, 0.025471556137328356, 0.00837663117351707, 0.0022749858369468675, 2.2566572299776048e-05, 1.6166079618194827e-07, 6.52263748461583e-09]
E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]


MC_flux = [0.05487097192089243, 0.04440512174328665, 0.018555746103953677, 0.019081170691700716, 0.04673392187154586, 0.021150826010103416, 0.035301110346496314, 0.15051425851121283, 0.3178388084844844, 0.2915480643163238]
MOC_flux = [0.055035098353831385, 0.04432925727674683, 0.017987406131318443, 0.017907561209162745, 0.04539360702480042, 0.019710012165971555, 0.03349040387916729, 0.1458221611019546, 0.3219260957418274, 0.29839839711521937]
MOC_tones_flux = [0.3113268304124537, 0.3367339423223383, 0.15368378967069207, 0.03589441503945606, 0.021196739568359668, 0.048475696277449326, 0.019170726124545286, 0.017976399353967117, 0.03651480575403351, 0.019026655476705016]
tones_flux = [6.52263748461583e-09, 1.6166079618194827e-07, 2.2566572299776048e-05, 0.0022749858369468675, 0.00837663117351707, 0.025471556137328356, 0.1040170872848463, 0.18740661552288226, 0.34047393864576864, 0.27220654333409494]

MC_total = [0.26481271, 0.41892899, 0.56278376, 0.68981225, 0.67524144, 0.48130316, 0.67820235, 0.85764397, 1.07211566, 1.76879026]
MC_absor = [0.01164354, 0.00897123, 0.07440102, 0.24198251, 0.28974487, 0.10185146, 0.29025338, 0.46792415, 0.68032592, 1.37187549]

tones_total  = [0.26205446, 0.42175232, 0.55982281, 0.65959835, 0.64546669, 0.47568846, 0.66879071, 0.86016966, 1.05421875, 5.04495411]
tones_absor  = [0.00956454, 0.00603311, 0.04377308, 0.14708432, 0.21239216, 0.02874255, 0.05339923, 0.09900321, 0.11985593, 0.77408611]
tones_nufis  = [0.03198346, 0.00929355, 0.0676909, 0.15872762, 0.11612148, 0.16457006, 0.55429328, 0.9051564, 1.32218379, 9.43462959]
tones_total.reverse()


MC_total.reverse()

diffEnergy = [E_bounds[i+1]-E_bounds[i] for i in range(10)]
MC_normalize = [MC_flux[i]*diffEnergy[i] for i in range(10)]
MOC_normalize = [MOC_flux[i]*diffEnergy[i] for i in range(10)]
MOC_tones_normalize = [MOC_tones_flux[i]*diffEnergy[i] for i in range(10)]
tones_normalize = [tones_flux[i]*diffEnergy[i] for i in range(10)]

print(MOC_tones_normalize)
print(tones_normalize)
MC_T_RR = []
MOC_T_RR = []
MOC_TONES_T_RR = []
TONES_T_RR = []
for i in range(10):
	MC_T_RR.append(MC_flux[i]*MC_total[i]*diffEnergy[i]/sum(MC_normalize))#[i])
	MOC_T_RR.append(MOC_flux[i]*MC_total[i]*diffEnergy[i]/sum(MC_normalize))#[i])
	MOC_TONES_T_RR.append(MOC_tones_flux[i]*tones_total[i]*diffEnergy[i]/sum(MOC_tones_normalize))#[i])
	TONES_T_RR.append(tones_flux[i]*tones_total[i]*diffEnergy[i]/sum(tones_normalize))#[i])

print(MOC_TONES_T_RR)
print(TONES_T_RR)
# print(MC_RR)

# plt.plot(E_bounds[:-1],MC_RR)

# plt.xscale('log')
# plt.show()
# E_bounds = [1e-5,0.058,0.14,0.28,0.625,4.,1e1,4e1,5.53e3,8.21e5,2.e7]
plt.step(E_bounds,[MC_T_RR[0]]+MC_T_RR,label='MC')
plt.step(E_bounds,[MOC_T_RR[0]]+MOC_T_RR,label='MOC using MC-generated XS')
plt.step(E_bounds,[MOC_TONES_T_RR[0]]+MOC_TONES_T_RR,label='MOC using Tones XS')
plt.step(E_bounds,[TONES_T_RR[0]]+TONES_T_RR,label='Tones Eq. Solution')
# plt.step(E_bounds,[tonesMOC[0]]+tonesMOC,label='MOC using Tones XS')
# plt.step(E_bounds,[tonesTheory[0]]+tonesTheory,label='Tones Eq. Solution')
plt.ylabel('Total Reaction Rate')
plt.xlabel('Energy (eV)')
plt.xscale('log')
plt.title('Total Reaction Rate in 6% enr. fuel pin using different methods')
plt.legend(loc='best')
plt.show()
