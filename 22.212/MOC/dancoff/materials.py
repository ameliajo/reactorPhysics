
class material:
    def __init__(self,SigmaT,SigmaF,SigmaS_matrix,chi):
        self.SigT = SigmaT
        self.SigF = SigmaF
        self.SigS_matrix = SigmaS_matrix
        self.chi  = chi

def getScatteringIntoG(gNow,SMatrix,phi):
    return sum([SMatrix[g_from][gNow]*phi[g_from] for g_from in range(len(SMatrix))])

def getFissionIntoG(gNow,chi,SigF,phi):
    return sum([chi[gNow]*SigF[g_from]*phi[g_from] for g_from in range(len(SigF))])









fuelTotal = [1.0e5]

modTotal = [1e-2]



fuelNuFission = [1.0e-1]
modNuFission = [0.0]

 # [[ a->a a->b a->c ],
 #  [ b->a b->b b->c ],
 #  [ c->a b->b c->c ]]

fuelScatter = [[0.0]]

modScatter  = [[1e-1]]

fuelChi = [1.0]
modChi = [0.0]





















