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
NVALS = 50000

minEnergy = 5.0
maxEnergy = 40.0

xvals = linspace(minEnergy, maxEnergy, num=NVALS)
E_width = xvals[1] - xvals[0]
capture, scatter = eval_xs(xvals, TEMP)


def plot_line(sigma):
    plt.plot([minEnergy,maxEnergy],[sigma,sigma])

def whichBandAmIIn(xsVal,bands):
    for i,band in enumerate(bands):
        if ((band[0] <= xsVal) and (xsVal <= band[1])):
            return i


def getProbabilities(xvals,xsVec,label):

    numBands = 5
    min_capture = np.amin(xsVec)-1e-3
    max_capture = np.amax(xsVec)+1e-3
    slice_spacing = (max_capture-min_capture)/float(numBands)
    bands = [(min_capture+i*slice_spacing,min_capture+(i+1)*slice_spacing) for i in range(0,numBands)]

    plot_line(bands[0][0])
    for band in bands:
        plot_line(band[1])

    currentBand = whichBandAmIIn(xsVec[0],bands)
    chunk_X = []
    chunk_Y = []

    probability = [0.0]*len(bands)
    
    for i in range(len(xvals)):
        if whichBandAmIIn(xsVec[i],bands) == currentBand:
            chunk_X.append(xvals[i])
            chunk_Y.append(xsVec[i])
        else:
            thisChunk = np.trapz(chunk_Y,chunk_X)/(chunk_X[-1]-chunk_X[0])
            probability[currentBand] = thisChunk
            currentBand = whichBandAmIIn(xsVec[i],bands)
            chunk_X = [xvals[i]]
            chunk_Y = [xsVec[i]]
    
    totalSum = sum(probability)
    probability = [val/totalSum for val in probability]
    plt.plot(xvals, xsVec, label=label)
    plt.show()
    return probability
print(getProbabilities(xvals,capture,'capture'))


