from pylab import *

def pdf(x):
	denom = 1 + (x - 3)**2
	return 1/pi/denom

def cdf(x):
	return 0.5 + arctan(x-3)/pi


xvals = linspace(-110, +110, 610)
figure()
plot(xvals, cdf(xvals), "g-", label="$CDF(x)$", alpha=0.75)
fill_between(xvals, zeros(len(xvals)), cdf(xvals), color='g', alpha=0.25)
plot(xvals, pdf(xvals), "k-", label="$PDF(x)$", linewidth=2)
xlim(-105, 105)
ylim(0, 1.1)
grid()
legend()
show()
