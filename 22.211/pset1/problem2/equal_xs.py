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

def plot_line(sigma):
    plt.plot([minEnergy,maxEnergy],[sigma,sigma])

def whichBandAmIIn(xsVal,bands):
    for i,band in enumerate(bands):
        if ((band[0] <= xsVal) and (xsVal <= band[1])):
            return i

def getProbabilities(xvals,xsVec,label,toPlot=True):
    numBands = 5
    min_capture = np.amin(xsVec)-1e-3
    max_capture = np.amax(xsVec)+1e-3
    slice_spacing = (max_capture-min_capture)/float(numBands)
    bands = [(min_capture+i*slice_spacing,min_capture+(i+1)*slice_spacing) for i in range(0,numBands)]


    currentBand = whichBandAmIIn(xsVec[0],bands)
    chunk_X = []
    chunk_Y = []

    xsIntegralPieces = [[] for i in range(len(bands))]

    if toPlot:
        plot_line(bands[0][0])
        for band in bands:
            plot_line(band[1])
        plt.plot(xvals, xsVec, label=label)
        plt.title('U-238 '+label+' cross section and bands')
        #plt.legend(loc='best')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Cross Section (b)')
        plt.yscale('log')
        plt.show()
    
    for i in range(len(xvals)):
        if whichBandAmIIn(xsVec[i],bands) == currentBand:
            chunk_X.append(xvals[i])
            chunk_Y.append(xsVec[i])
        else:
            # We are moving into a new band, so want to integrate the progress
            # we made thus far and add it to our tally vector
            thisChunk = np.trapz(chunk_Y,chunk_X)#/(chunk_X[-1]-chunk_X[0])
            xsIntegralPieces[currentBand].append(thisChunk)
            currentBand = whichBandAmIIn(xsVec[i],bands)
            chunk_X = [xvals[i]]
            chunk_Y = [xsVec[i]]
    
    return xsIntegralPieces

TEMP = 0
NVALS = 50000

minEnergy = 5.0
maxEnergy = 40.0

xvals = linspace(minEnergy, maxEnergy, num=NVALS)
capture, scatter = eval_xs(xvals, TEMP)
total = [capture[i]+scatter[i] for i in range(len(scatter))]

xsIntegralPieces = getProbabilities(xvals,total,'total',True)
probabilities = [0.0]*len(xsIntegralPieces)
for i,val in enumerate(xsIntegralPieces):
    probabilities[i] = sum(xsIntegralPieces[i])/len(xsIntegralPieces[i])

sumVal = sum(probabilities)
probabilities = [val/sumVal for val in probabilities]
print(probabilities)




