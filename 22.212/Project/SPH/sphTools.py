

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
        #sph_mod  = MC_modFlux[0][g] / MOC_modFlux[0][g]
        #sph_fuel = MC_fuelFlux[0][g] / MOC_fuelFlux[0][g]
        sph_mod  = MOC_modFlux[0][g] / MC_modFlux[0][g]
        sph_fuel = MOC_fuelFlux[0][g] / MC_fuelFlux[0][g]

        sph[0].append(sph_mod)
        sph[1].append(sph_fuel)


    return sph







