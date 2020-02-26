import matplotlib.pyplot as plt
import numpy as np
from math import pi, exp

def getAlphaRange(beta,A,kbT,E,nAlpha):
    aMin = ( (E)**0.5 - (E+beta*kbT)**0.5 )**2 / (A*kbT)
    aMax = ( (E)**0.5 + (E+beta*kbT)**0.5 )**2 / (A*kbT)
    return np.linspace(aMin,aMax,nAlpha)

def getFreeGasSAB(alpha,beta):
    return (4.0*pi*alpha)**-0.5 * exp( -(alpha+beta)**2 / (4.0*alpha) )

E    = 1.0
A    = 1.0
kbT  = 8.617e-5 * 296.0
beta = 10.0

def getAlphaPDF(beta,kbT,A,E,numAlpha=50000):
    alphaRange = getAlphaRange(beta,A,kbT,E,numAlpha)
    sabChunk = [ getFreeGasSAB(alpha,beta) for alpha in alphaRange ]
    scaling = 1.0/np.trapz(sabChunk,alphaRange)
    eq15 = [sabChunk[a]*scaling for a in range(len(sabChunk))]
    return alphaRange, eq15

def getAlphaCDF(alphaRange, alphaPDF):
    alphaCDF = [0.0]*len(alphaRange)
    for a in range(1,len(alphaRange)):
        alphaCDF[a] = alphaCDF[a-1]+alphaPDF[a]
    invTotal = 1.0/alphaCDF[len(alphaCDF)-1]
    return [CDF*invTotal for CDF in alphaCDF]

alphaRange, alphaPDF = getAlphaPDF(beta,kbT,A,E)
plt.xlabel('alpha')
plt.title('H in H2O Free Gas S(a,b) alpha PDF (296 K, beta=10)'); 
plt.plot(alphaRange,alphaPDF,label='alpha PDF'); plt.show()
alphaCDF = getAlphaCDF(alphaRange,alphaPDF)
plt.xlabel('alpha')
plt.title('H in H2O Free Gas S(a,b) alpha CDF (296 K, beta=10)'); 
plt.plot(alphaRange,alphaCDF,label='alpha CDF'); plt.show()


 



