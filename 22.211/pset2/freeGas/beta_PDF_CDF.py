import matplotlib.pyplot as plt
import numpy as np
from math import pi, exp




def getAlphaRange(beta,A,kbT,E,nAlpha):
    aMin = ( (E)**0.5 - (E+beta*kbT)**0.5 )**2 / (A*kbT)
    aMax = ( (E)**0.5 + (E+beta*kbT)**0.5 )**2 / (A*kbT)
    return np.linspace(aMin,aMax,nAlpha)


def getFreeGasSAB(alpha,beta):
    return (4.0*pi*alpha)**-0.5 * exp( -(alpha+beta)**2 / (4.0*alpha) )




def getAlphaPDF(beta,kbT,A,E):
    alphaRange = getAlphaRange(beta,A,kbT,E,50000)
    sabChunk = [ getFreeGasSAB(alpha,beta) for alpha in alphaRange ]
    denominator = np.trapz(sabChunk,alphaRange)
    eq15 = [sabChunk[a]/denominator for a in range(len(sabChunk))]
    print(np.trapz(eq15,alphaRange))
    return alphaRange, eq15

def getAlphaCDF(alphaRange, alphaPDF):
    alphaCDF = [0.0]*len(alphaRange)
    for a in range(1,len(alphaRange)):
        alphaCDF[a] = alphaCDF[a-1]+alphaPDF[a]
    invTotal = 1.0/alphaCDF[len(alphaCDF)-1]
    return [CDF*invTotal for CDF in alphaCDF]



E = 1.0
kb = 8.617e-5
kbT = kb * 296.0
A = 1.0

betaMin = -E/kbT
betaMax = 20.0

betas = list(np.linspace(betaMin,betaMax,1000)) 
betas.sort()

eq14 = [0.0]*len(betas)
for b,beta in enumerate(betas):
    alphaRange = getAlphaRange(beta,A,kbT,E,2000)
    sabChunk = [ getFreeGasSAB(alpha,beta) for alpha in alphaRange ]
    eq14[b] = np.trapz(sabChunk,alphaRange)

    
eq14Denominator = np.trapz(eq14,betas)
eq14 = [PDFVal/eq14Denominator for PDFVal in eq14]

cdf = [0.0]*len(eq14)
for b in range(1,len(eq14)):
    cdf[b] = cdf[b-1]+eq14[b]
cdf = [x/cdf[len(cdf)-1] for x in cdf]
plt.plot(betas,eq14)
plt.plot(betas,cdf)
plt.show()

