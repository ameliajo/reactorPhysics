import matplotlib.pyplot as plt
import numpy as np
from math import pi, exp

def getAlphaRange(beta,A,kbT,E,nAlpha):
    aMin = ( (E)**0.5 - (E+beta*kbT)**0.5 )**2 / (A*kbT)
    aMax = ( (E)**0.5 + (E+beta*kbT)**0.5 )**2 / (A*kbT)
    return np.linspace(aMin,aMax,nAlpha)

def getFreeGasSAB(alpha,beta):
    return (4.0*pi*alpha)**-0.5 * exp( -(alpha+beta)**2 / (4.0*alpha) )

def getBetaPDF(A,kbT,E):
    betaMin, betaMax = -E/kbT, 20.0
    betas = list(np.linspace(betaMin,betaMax,1000)) 
    eq14 = [0.0]*len(betas)
    for b,beta in enumerate(betas):
        alphaRange = getAlphaRange(beta,A,kbT,E,2000)
        sabChunk = [ getFreeGasSAB(alpha,beta) for alpha in alphaRange ]
        eq14[b] = np.trapz(sabChunk,alphaRange)
    eq14Denominator = np.trapz(eq14,betas)
    eq14 = [PDFVal/eq14Denominator for PDFVal in eq14]
    return betas, eq14

def getBetaCDF(betas,betaPDF):
    betaCDF = [0.0]*len(betaPDF)
    for b in range(1,len(betaPDF)):
        betaCDF[b] = betaCDF[b-1]+betaPDF[b]
    return [x/betaCDF[-1] for x in betaCDF]

E = 1.0
A = 1.0
kbT = 8.617e-5 * 296.0

betas, betaPDF = getBetaPDF(A,kbT,E)
plt.title('H in H2O Free Gas S(a,b) beta PDF (296 K)'); 
plt.plot(betas,betaPDF)
plt.show()

betaCDF = getBetaCDF(betas,betaPDF)
plt.title('H in H2O Free Gas S(a,b) beta CDF (296 K)'); 
plt.plot(betas,betaCDF)
plt.show()

