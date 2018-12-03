

def normalizeToUnity(nPins, modFlux, fuelFlux):
    #invSum = 1.0/sum( [sum(m) for m in modFlux] + [sum(f) for f in fuelFlux] )
    #for i in range(nPins):
    #    modFlux[i]  = [ x * invSum for x in modFlux[i]  ]
    #    fuelFlux[i] = [ x * invSum for x in fuelFlux[i] ]

    invSum = 1.0/sum( [sum(m) for m in modFlux] )
    for i in range(nPins): modFlux[i]  = [ x * invSum for x in modFlux[i]  ]

    invSum = 1.0/sum( [sum(m) for m in fuelFlux] )
    for i in range(nPins): fuelFlux[i]  = [ x * invSum for x in fuelFlux[i]  ]



def calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux):
    sph = [[],[]]
    for g in range(nGroups):
        sph[0].append(MOC_modFlux[0][g] / MC_modFlux[0][g])
        sph[1].append(MOC_fuelFlux[0][g] / MC_fuelFlux[0][g])

    return sph


def readInEigenvalue(filename):
    for line in reversed(list(open(filename))):
        if 'Combined k-effective' in line:
            return line.split()[3]



def plotThis(MC_modFlux,MC_fuelFlux):
    plt.plot(MC_modFlux[0],label='MC mod')
    plt.plot(MOC_modFlux[0],label='MOC mod')
    plt.plot(MC_fuelFlux[0],label='MC fuel')
    plt.plot(MOC_fuelFlux[0],label='MOC fuel')
    plt.legend(loc='best')




