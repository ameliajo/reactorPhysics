# 22.211/PSet01/Problem1
#
# Doppler Broadening of SLBW resonance

from resonance import all_faddeeva_resonances as RESONANCES
from pylab import *

TEMP = 1000
TEMP = 0
NVALS = 10000
E0 = 1E-2
E1 = 1E+2

xvals = logspace(-2, +2, num=NVALS)
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

# Plotting parameters
fig, [ax1, ax2] = subplots(2, sharex=True, sharey=True)
ax1.loglog(xvals, y_capture, 'b-', label="capture")
ax2.loglog(xvals, y_scatter, 'r-', label="scatter")
ax1.set_ylabel("$\sigma_\gamma$ (barns)", fontsize=12)
ax2.set_ylabel("$\sigma_n$ (barns)", fontsize=12)
ax1.set_title("Problem 1: U-238 SLBW Cross Sections at {} K".format(TEMP), fontweight="bold")
ax2.set_xlabel("Energy (eV)", fontsize=12)
tight_layout()
figure(2)
loglog(xvals, y_capture, 'b-', label="capture")
loglog(xvals, y_scatter, 'r-', label="scatter")
loglog(xvals, y_capture+y_scatter, '-', c='purple', label="total")
xlabel("Energy (eV)", fontsize=12)
ylabel("$\sigma$ (barns)", fontsize=12)
title("Problem 1: U-238 SLBW Cross Section at {} K".format(TEMP), fontweight="bold")
legend()
show()
