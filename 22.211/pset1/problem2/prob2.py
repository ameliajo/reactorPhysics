# 22.211/PSet01/Problem1
#
# Doppler Broadening of SLBW resonance

from resonance import all_faddeeva_resonances as RESONANCES
from pylab import *
import numpy as np

def eval_xs(energies, TEMP):
    y_capture = zeros(NVALS, dtype=float)
    y_scatter = zeros(NVALS, dtype=float)
    for i, e in enumerate(xvals):
    	for res in RESONANCES:
    		sigma_y = res.sigma_y(e, t=TEMP)
    		y_capture[i] += sigma_y
    		sigma_n = res.sigma_n(e, t=TEMP)
    		y_scatter[i] += sigma_n
    # Don't count the potential xs multiple time
    y_scatter -= res.potential*(len(RESONANCES)-1)
    return y_capture, y_scatter

TEMP = 0
NVALS = 10000
E0 = 1E-2
E1 = 1E+2
xvals = linspace(5, 40, num=NVALS)
E_width = xvals[1] - xvals[0]
capture, scatter = eval_xs(xvals, TEMP)

# 0 K Evaluation of multigroup cross sections
capture_group_xs_0K = np.zeros([5,])
scatter_group_xs_0K = np.zeros([5,])

# Find bounds for problem
capture_max = np.amax(capture)
capture_bounds = np.linspace(0, capture_max, 6)
scatter_max = np.amax(scatter)
scatter_bounds = np.linspace(0, scatter_max, 6)
for i, e in enumerate(xvals):
    cap_i = capture[i]
    c_idx = np.searchsorted(capture_bounds, cap_i) - 1
    capture_group_xs_0K[c_idx] += E_width*cap_i
    scat_i = scatter[i]
    s_idx = np.searchsorted(scatter_bounds, scat_i) - 1
    scatter_group_xs_0K[c_idx] += E_width*scat_i

print('T = {:f} K \nE_width = {:f} \nBounds(barns): Capture, Bounds(barns): Scatter'.format(TEMP, E_width))
for i in range(5):
    print('{:08.2f}-{:08.2f}: {:08.2f}, {:08.2f}-{:08.2f}: {:08.2f},'.format(\
          capture_bounds[i], capture_bounds[i+1], capture_group_xs_0K[i],
          scatter_bounds[i], scatter_bounds[i+1], scatter_group_xs_0K[i]))

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.loglog(xvals, capture, label='capture')
ax1.loglog(xvals, scatter, label='scatter')
ax1.set_xlabel('Energy (eV)')
ax1.set_ylabel('Barns')
for i in range(5):
    ax1.axhline(y=capture_bounds[i+1], color='r', label='Capture Bound')
    ax1.axhline(y=scatter_bounds[i+1], color='g', label='Scatter Bound')
ax1.legend()
ax1.set_title('0 K')

TEMP = 1000
NVALS = 10000
E0 = 1E-2
E1 = 1E+2
xvals = linspace(5, 40, num=NVALS)
E_width = xvals[1] - xvals[0]
capture, scatter = eval_xs(xvals, TEMP)

# 0 K Evaluation of multigroup cross sections
capture_group_xs_0K = np.zeros([5,])
scatter_group_xs_0K = np.zeros([5,])

# Find bounds for problem
capture_max = np.amax(capture)
capture_bounds = np.linspace(0, capture_max, 6)
scatter_max = np.amax(scatter)
scatter_bounds = np.linspace(0, scatter_max, 6)
for i, e in enumerate(xvals):
    cap_i = capture[i]
    c_idx = np.searchsorted(capture_bounds, cap_i) - 1
    capture_group_xs_0K[c_idx] += E_width*cap_i
    scat_i = scatter[i]
    s_idx = np.searchsorted(scatter_bounds, scat_i) - 1
    scatter_group_xs_0K[c_idx] += E_width*scat_i

print('T = {:f} K \nE_width = {:f} \nBounds(barns): Capture, Bounds(barns): Scatter'.format(TEMP, E_width))
for i in range(5):
    print('{:08.2f}-{:08.2f}: {:08.2f}, {:08.2f}-{:08.2f}: {:08.2f},'.format(\
          capture_bounds[i], capture_bounds[i+1], capture_group_xs_0K[i],
          scatter_bounds[i], scatter_bounds[i+1], scatter_group_xs_0K[i]))

ax2 = fig.add_subplot(212)
ax2.loglog(xvals, capture, label='capture')
ax2.loglog(xvals, scatter, label='scatter')
ax2.set_xlabel('Energy (eV)')
ax2.set_ylabel('Barns')
for i in range(5):
    ax2.axhline(y=capture_bounds[i+1], color='r', label='Capture Bound')
    ax2.axhline(y=scatter_bounds[i+1], color='g', label='Scatter Bound')
ax2.legend()
ax2.set_title('1000 K')

plt.show()




# Plotting parameters
# fig, [ax1, ax2] = subplots(2, sharex=True, sharey=True)
# ax1.loglog(xvals, y_capture, 'b-', label="capture")
# ax2.loglog(xvals, y_scatter, 'r-', label="scatter")
# ax1.set_ylabel("$\sigma_\gamma$ (barns)", fontsize=12)
# ax2.set_ylabel("$\sigma_n$ (barns)", fontsize=12)
# ax1.set_title("Problem 1: Cross Sections at {} K".format(TEMP), fontweight="bold")
# ax2.set_xlabel("Energy (eV)", fontsize=12)
# tight_layout()
# figure(2)
# loglog(xvals, y_capture, 'b-', label="capture")
# loglog(xvals, y_scatter, 'r-', label="scatter")
# loglog(xvals, y_capture+y_scatter, '-', c='purple', label="total")
# xlabel("Energy (eV)", fontsize=12)
# ylabel("$\sigma$ (barns)", fontsize=12)
# title("Problem 1: Cross Section at {} K".format(TEMP), fontweight="bold")
# legend()
# show()
