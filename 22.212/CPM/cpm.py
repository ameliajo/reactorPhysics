from numpy import array
from numpy.linalg import norm
import matplotlib.pyplot as plt
from data import *
from fuel_matrix_data import *
from mod_matrix_data import *
from pvvc_data import *

scat = [ scat_f0, scat_f1, scat_f2, scat_f3, scat_f4, scat_f5, scat_f6, scat_f7, scat_f8, 
         scat_m0, scat_m1, scat_m2, scat_m3, scat_m4, scat_m5, scat_m6, scat_m7, scat_m8 ]

num_groups = 10
num_mat = 18
num_iter = 0
# set phi_0 and k_0

# phi_f1 = [ fuel_1_group_1 fuel_1_group_2 fuel_1_group_3 ... ]
# phi_f2 = [ fuel_2_group_1 fuel_2_group_2 fuel_2_group_3 ... ]
#         ...
# phi_m1 = [ mod__1_group_1 mod__1_group_2 mod__1_group_3 ... ]
# phi_m2 = [ mod__2_group_1 mod__2_group_2 mod__2_group_3 ... ]

phi_old = []
phi_new = []
for entry in range(18):
	phi_old.append([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
	phi_new.append([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])

# Initialize to dummy values
k = 1
Q = [0]*num_mat
converged = False

while (not converged):
	for g in range(num_groups):
		# Q = Sigma_s * phi + nu Sigma_f * phi / k
		# where sigma_s is sum over g' of all g'->g * phi_g'
		#      scat 0 -> 0 * phi 0 + scat 1 -> 0 * phi 1

		Q = [0] * num_mat
		for g_prime in range(num_groups):
			for mat in range(num_mat):
				#         S g' -> g * phi g'                      + Chi(g) * nu * Sigma_f g' * phi g'
				Q[mat] += scat[mat][g_prime][g]*phi_old[mat][g_prime] + chi[mat][g]*nuf[mat][g_prime]*phi_old[mat][g_prime]

		#  phi_new = PvvC * Q
		#          = all the probability of getting to this cell * getting in this energy group
		# pvvc_E0[a][b] is the probability of going from mat a --> mat b in energy group 0

		# so whats the probability of getting from x --> f0, for all x?

		pvvc = [0]*num_mat

		pvvc_Eg = pvvc_all[g]

		for mat1 in range(num_mat):
			for mat2 in range(num_mat):
				phi_new[mat1][g] += pvvc_Eg[mat2][mat1]*Q[mat2]
				# Trying to find all the ways that a neutron could go from material 2 --> material 1


	# k_new = Sum of nu sigma f phi_1 / sum of nu sigma f phi_0
	k_new_numerator = 0
	k_new_denominator = 0
	for energy in range(num_groups):
		for mat in range(num_mat):
			k_new_numerator   += nuf[mat][energy]*phi_new[mat][energy]
			k_new_denominator += nuf[mat][energy]*phi_old[mat][energy]
	k_new = k_new_numerator / k_new_denominator

	phi_diff = 0.0

	for i in range(len(phi_old)):
		for j in range(len(phi_old[i])):
			phi_old[i][j] /= norm(phi_old[i])
			phi_new[i][j] /= norm(phi_new[i])
			phi_diff = abs(phi_old[i][j]-phi_new[i][j])

	num_iter += 1
	if abs(k_new-k)/k < 1e-5 or phi_diff < 1e-1:
		converged = True

	k = k_new

	for i in range(len(phi_old)):
		phi_old[i] = phi_new[i][:]


	# for i in range(len(phi_old)):
	# 	for j in range(len(phi_old[i])):
	# 		phi_old[i][j] /= norm(phi_old[i])


print("k ",k,"num iterations",num_iter)
# print(phi_new[0])
for i in range(9):
	plt.plot(phi_old[i])
plt.show()














# # Set geometry and material properties


# for g in range(num_groups):
# 	# Load collision probabilities

# 	# Normalize CPs using sybrhl.m (transpose your Σt vector)

# 	# Extract individual matrices Pvv , Pvs , Psv , Pss

# 	# Build the PCvv matrix with boundary matrix

# 	# Guess initial k and flux, and calculate sources

# 	# Calculate flux and k

# 	# Iterate to convergence on both flux and k
