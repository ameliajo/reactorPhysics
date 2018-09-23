from data import *
from fuel_matrix_data import *
from mod_matrix_data import *
from pvvc_data import *


num_groups = 10
num_mat = 18
# set phi_0 and k_0

# phi_f1 = [ fuel_1_group_1 fuel_1_group_2 fuel_1_group_3 ... ]
# phi_f2 = [ fuel_2_group_1 fuel_2_group_2 fuel_2_group_3 ... ]
#         ...
# phi_m1 = [ mod__1_group_1 mod__1_group_2 mod__1_group_3 ... ]
# phi_m2 = [ mod__2_group_1 mod__2_group_2 mod__2_group_3 ... ]
# 

phi_f0 = [1]*num_groups
phi_f1 = [1]*num_groups
phi_f2 = [1]*num_groups
phi_f3 = [1]*num_groups
phi_f4 = [1]*num_groups
phi_f5 = [1]*num_groups
phi_f6 = [1]*num_groups
phi_f7 = [1]*num_groups
phi_f8 = [1]*num_groups

phi_m0 = [1]*num_groups
phi_m1 = [1]*num_groups
phi_m2 = [1]*num_groups
phi_m3 = [1]*num_groups
phi_m4 = [1]*num_groups
phi_m5 = [1]*num_groups
phi_m6 = [1]*num_groups
phi_m7 = [1]*num_groups
phi_m8 = [1]*num_groups


phi_f0_new = [1]*num_groups
phi_f1_new = [1]*num_groups
phi_f2_new = [1]*num_groups
phi_f3_new = [1]*num_groups
phi_f4_new = [1]*num_groups
phi_f5_new = [1]*num_groups
phi_f6_new = [1]*num_groups
phi_f7_new = [1]*num_groups
phi_f8_new = [1]*num_groups

phi_m0_new = [1]*num_groups
phi_m1_new = [1]*num_groups
phi_m2_new = [1]*num_groups
phi_m3_new = [1]*num_groups
phi_m4_new = [1]*num_groups
phi_m5_new = [1]*num_groups
phi_m6_new = [1]*num_groups
phi_m7_new = [1]*num_groups
phi_m8_new = [1]*num_groups


# Initialize to dummy values
k = 1
Q   = [0]*num_mat
converged = False

while (not converged):
	for g in range(num_groups):
		# Q = Sigma_s * phi + nu Sigma_f * phi / k
		# where sigma_s is sum over g' of all g'->g * phi_g'
		#      scat 0 -> 0 * phi 0 + scat 1 -> 0 * phi 1
		for g_prime in range(num_groups):
			# Q[mat] += S g' -> g * phi g'              + Chi(g) * nu*Sigma_f g' * phi g'
			Q[0] += scat_f0[g_prime][g]*phi_f0[g_prime] + chi_f0[g]*nuf_f0[g_prime]*phi_f0[g_prime]
			Q[1] += scat_f1[g_prime][g]*phi_f1[g_prime] + chi_f1[g]*nuf_f1[g_prime]*phi_f1[g_prime]
			Q[2] += scat_f2[g_prime][g]*phi_f2[g_prime] + chi_f2[g]*nuf_f2[g_prime]*phi_f2[g_prime]
			Q[3] += scat_f3[g_prime][g]*phi_f3[g_prime] + chi_f3[g]*nuf_f3[g_prime]*phi_f3[g_prime]
			Q[4] += scat_f4[g_prime][g]*phi_f4[g_prime] + chi_f4[g]*nuf_f4[g_prime]*phi_f4[g_prime]
			Q[5] += scat_f5[g_prime][g]*phi_f5[g_prime] + chi_f5[g]*nuf_f5[g_prime]*phi_f5[g_prime]
			Q[6] += scat_f6[g_prime][g]*phi_f6[g_prime] + chi_f6[g]*nuf_f6[g_prime]*phi_f6[g_prime]
			Q[7] += scat_f7[g_prime][g]*phi_f7[g_prime] + chi_f7[g]*nuf_f7[g_prime]*phi_f7[g_prime]
			Q[8] += scat_f8[g_prime][g]*phi_f8[g_prime] + chi_f8[g]*nuf_f8[g_prime]*phi_f8[g_prime]

			Q[9]  += scat_m0[g_prime][g]*phi_m0[g_prime] + chi_m0[g]*nuf_m0[g_prime]*phi_m0[g_prime]
			Q[10] += scat_m1[g_prime][g]*phi_m1[g_prime] + chi_m1[g]*nuf_m1[g_prime]*phi_m1[g_prime]
			Q[11] += scat_m2[g_prime][g]*phi_m2[g_prime] + chi_m2[g]*nuf_m2[g_prime]*phi_m2[g_prime]
			Q[12] += scat_m3[g_prime][g]*phi_m3[g_prime] + chi_m3[g]*nuf_m3[g_prime]*phi_m3[g_prime]
			Q[13] += scat_m4[g_prime][g]*phi_m4[g_prime] + chi_m4[g]*nuf_m4[g_prime]*phi_m4[g_prime]
			Q[14] += scat_m5[g_prime][g]*phi_m5[g_prime] + chi_m5[g]*nuf_m5[g_prime]*phi_m5[g_prime]
			Q[15] += scat_m6[g_prime][g]*phi_m6[g_prime] + chi_m6[g]*nuf_m6[g_prime]*phi_m6[g_prime]
			Q[16] += scat_m7[g_prime][g]*phi_m7[g_prime] + chi_m7[g]*nuf_m7[g_prime]*phi_m7[g_prime]
			Q[17] += scat_m8[g_prime][g]*phi_m8[g_prime] + chi_m8[g]*nuf_m8[g_prime]*phi_m8[g_prime]

		#  phi_new = PvvC * Q
		#          = all the probability of getting to this cell * getting in this energy group
		# pvvc_E0[a][b] is the probability of going from mat a --> mat b in energy group 0

		# so whats the probability of getting from x --> f0, for all x?
		pvvc_f0_Eg = 0
		pvvc_f1_Eg = 0
		pvvc_f2_Eg = 0
		pvvc_f3_Eg = 0
		pvvc_f4_Eg = 0
		pvvc_f5_Eg = 0
		pvvc_f6_Eg = 0
		pvvc_f7_Eg = 0
		pvvc_f8_Eg = 0
		
		pvvc_m0_Eg = 0
		pvvc_m1_Eg = 0
		pvvc_m2_Eg = 0
		pvvc_m3_Eg = 0
		pvvc_m4_Eg = 0
		pvvc_m5_Eg = 0
		pvvc_m6_Eg = 0
		pvvc_m7_Eg = 0
		pvvc_m8_Eg = 0

		pvvc_Eg = pvvc_all[g]

		for material in range(len(pvvc_Eg)):
			pvvc_f0_Eg += pvvc_Eg[material][0]
			pvvc_f1_Eg += pvvc_Eg[material][1]
			pvvc_f2_Eg += pvvc_Eg[material][2]
			pvvc_f3_Eg += pvvc_Eg[material][3]
			pvvc_f4_Eg += pvvc_Eg[material][4]
			pvvc_f5_Eg += pvvc_Eg[material][5]
			pvvc_f6_Eg += pvvc_Eg[material][6]
			pvvc_f7_Eg += pvvc_Eg[material][7]
			pvvc_f8_Eg += pvvc_Eg[material][8]

			pvvc_m0_Eg += pvvc_Eg[material][9]
			pvvc_m1_Eg += pvvc_Eg[material][10]
			pvvc_m2_Eg += pvvc_Eg[material][11]
			pvvc_m3_Eg += pvvc_Eg[material][12]
			pvvc_m4_Eg += pvvc_Eg[material][13]
			pvvc_m5_Eg += pvvc_Eg[material][14]
			pvvc_m6_Eg += pvvc_Eg[material][15]
			pvvc_m7_Eg += pvvc_Eg[material][16]
			pvvc_m8_Eg += pvvc_Eg[material][17]



		phi_f0_new[g] = pvvc_f0_Eg * Q[0]
		phi_f1_new[g] = pvvc_f1_Eg * Q[1]
		phi_f2_new[g] = pvvc_f2_Eg * Q[2]
		phi_f3_new[g] = pvvc_f3_Eg * Q[3]
		phi_f4_new[g] = pvvc_f4_Eg * Q[4]
		phi_f5_new[g] = pvvc_f5_Eg * Q[5]
		phi_f6_new[g] = pvvc_f6_Eg * Q[6]
		phi_f7_new[g] = pvvc_f7_Eg * Q[7]
		phi_f8_new[g] = pvvc_f8_Eg * Q[8]

		phi_m0_new[g] = pvvc_m0_Eg * Q[9]
		phi_m1_new[g] = pvvc_m1_Eg * Q[10]
		phi_m2_new[g] = pvvc_m2_Eg * Q[11]
		phi_m3_new[g] = pvvc_m3_Eg * Q[12]
		phi_m4_new[g] = pvvc_m4_Eg * Q[13]
		phi_m5_new[g] = pvvc_m5_Eg * Q[14]
		phi_m6_new[g] = pvvc_m6_Eg * Q[15]
		phi_m7_new[g] = pvvc_m7_Eg * Q[16]
		phi_m8_new[g] = pvvc_m8_Eg * Q[17]

	# k_new = Sum of nu sigma f phi_1 * k0 / sum of nu sigma f phi_0
	print(phi_f0_new)
	k_new_numerator = 0
	k_new_denominator = 0
	for energy in range(num_groups):
		k_new_numerator   += nuf_f0[energy]*phi_f0_new[energy]
		k_new_denominator += nuf_f0[energy]*phi_f0[energy]
	k_new = k_new_numerator * k / k_new_denominator
	print(k)
	
	if abs(k_new-k)/k < 1e-5:
		converged = True

	k = k_new

	
	phi_f0 = phi_f0_new	
	phi_f1 = phi_f1_new	
	phi_f2 = phi_f2_new	
	phi_f3 = phi_f3_new	
	phi_f4 = phi_f4_new	
	phi_f5 = phi_f5_new	
	phi_f6 = phi_f6_new	
	phi_f7 = phi_f7_new	
	phi_f8 = phi_f8_new	

	phi_m0 = phi_m0_new	
	phi_m1 = phi_m1_new	
	phi_m2 = phi_m2_new	
	phi_m3 = phi_m3_new	
	phi_m4 = phi_m4_new	
	phi_m5 = phi_m5_new	
	phi_m6 = phi_m6_new	
	phi_m7 = phi_m7_new	
	phi_m8 = phi_m8_new	



















# # Set geometry and material properties


# for g in range(num_groups):
# 	# Load collision probabilities

# 	# Normalize CPs using sybrhl.m (transpose your Σt vector)

# 	# Extract individual matrices Pvv , Pvs , Psv , Pss

# 	# Build the PCvv matrix with boundary matrix

# 	# Guess initial k and flux, and calculate sources

# 	# Calculate flux and k

# 	# Iterate to convergence on both flux and k
