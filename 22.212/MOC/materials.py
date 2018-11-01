
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

