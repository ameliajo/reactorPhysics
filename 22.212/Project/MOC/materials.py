
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









fuelTotal  = [2.73e-1, 4.31e-1]
fuelTotal1 = [2.73e-1, 4.31e-1]
fuelTotal2 = [2.73e-1, 4.31e-1]
fuelTotal3 = [2.73e-1, 4.31e-1]
fuelTotal4 = [2.73e-1, 4.31e-1]
fuelTotal5 = [2.73e-1, 4.31e-1]
fuelTotal6 = [2.73e-1, 4.31e-1]
fuelTotal7 = [2.73e-1, 4.31e-1]
fuelTotal8 = [2.73e-1, 4.31e-1]
fuelTotal9 = [2.73e-1, 4.31e-1]

modTotal   = [2.01e-1, 6.53e-1]
modTotal1  = [2.01e-1, 6.53e-1]
modTotal2  = [2.01e-1, 6.53e-1]
modTotal3  = [2.01e-1, 6.53e-1]
modTotal4  = [2.01e-1, 6.53e-1]
modTotal5  = [2.01e-1, 6.53e-1]
modTotal6  = [2.01e-1, 6.53e-1]
modTotal7  = [2.01e-1, 6.53e-1]
modTotal8  = [2.01e-1, 6.53e-1]
modTotal9  = [2.01e-1, 6.53e-1]



fuelNuFission  = [2.58e-2, 1.51e-3]
fuelNuFission1 = [2.58e-2, 1.51e-3]
fuelNuFission2 = [2.58e-2, 1.51e-3]
fuelNuFission3 = [2.58e-2, 1.51e-3]
fuelNuFission4 = [2.58e-2, 1.51e-3]
fuelNuFission5 = [2.58e-2, 1.51e-3]
fuelNuFission6 = [2.58e-2, 1.51e-3]
fuelNuFission7 = [2.58e-2, 1.51e-3]
fuelNuFission8 = [2.58e-2, 1.51e-3]
fuelNuFission9 = [2.58e-2, 1.51e-3]

modNuFission  = [0.00000, 0.00000]
modNuFission1  = [0.00000, 0.00000]
modNuFission2  = [0.00000, 0.00000]
modNuFission3  = [0.00000, 0.00000]
modNuFission4  = [0.00000, 0.00000]
modNuFission5  = [0.00000, 0.00000]
modNuFission6  = [0.00000, 0.00000]
modNuFission7  = [0.00000, 0.00000]
modNuFission8  = [0.00000, 0.00000]
modNuFission9  = [0.00000, 0.00000]

 # [[ a->a a->b a->c ],
 #  [ b->a b->b b->c ],
 #  [ c->a b->b c->c ]]

fuelScatter = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter1 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter2 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter3 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter4 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter5 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter6 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter7 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter8 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]
fuelScatter9 = [[1.33e-01,4.95e-02],
               [0.00e+00,3.79e-01,]]






modScatter  = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter1 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter2 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter3 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter4 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter5 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter6 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter7 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter8 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]
modScatter9 = [[4.97e-02,8.20e-02],
               [0.00e+00,1.84e-01]]












fuelChi = [7.0e-1, 3.0e-1]
fuelChi1= [7.0e-1, 3.0e-1]
fuelChi2= [7.0e-1, 3.0e-1]
fuelChi3= [7.0e-1, 3.0e-1]
fuelChi4= [7.0e-1, 3.0e-1]
fuelChi5= [7.0e-1, 3.0e-1]
fuelChi6= [7.0e-1, 3.0e-1]
fuelChi7= [7.0e-1, 3.0e-1]
fuelChi8= [7.0e-1, 3.0e-1]
fuelChi9= [7.0e-1, 3.0e-1]

modChi1 = [0.0000, 0.0000]
modChi2 = [0.0000, 0.0000]
modChi3 = [0.0000, 0.0000]
modChi4 = [0.0000, 0.0000]
modChi5 = [0.0000, 0.0000]
modChi6 = [0.0000, 0.0000]
modChi7 = [0.0000, 0.0000]
modChi8 = [0.0000, 0.0000]
modChi9 = [0.0000, 0.0000]





















