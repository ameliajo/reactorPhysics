import numpy as np




def getQ(cell):
    return cell.src + 0.5*cell.SigS*cell.phi


def step(psi_in,mu_n,cell,dx):
    Q = getQ(cell)
    mu = abs(mu_n)
    return (dx*Q/mu) + psi_in * (1.0 - dx*cell.SigT/mu)


def diamondDiff(psi_in,mu_n,cell,dx):
    Q = getQ(cell)
    mu = abs(mu_n)
    return ( ( 2.0*mu - cell.SigT*dx ) * psi_in + 2.0*Q*dx ) / \
             ( 2.0*mu + cell.SigT*dx )

def stepCharacteristic(psi_in,mu_n,cell,dx):
    Q = getQ(cell)
    mu = abs(mu_n)
    return (1.0 - (dx*cell.SigT/mu)*np.exp(-cell.SigT*dx/mu)) * psi_in + \
           ((dx*Q/mu)*np.exp(-cell.SigT*dx/mu))


def getPsiOut(psi_in,mu,cell,dx,method):
    if method == 'step':               return step(psi_in,mu,cell,dx)
    if method == 'diamond':            return diamondDiff(psi_in,mu,cell,dx)
    if method == 'stepCharacteristic': return stepCharacteristic(psi_in,mu,cell,dx)
    raise ValueError

        


def broom(S,cells,dx,method):
    assert(S.N%2 == 0) # Require that order is even

    psi = np.zeros((len(cells)+1,S.N))

    # Forward Sweep
    for n in range(int(S.N/2)):
        mu = S.mu[n]
        psi_in = 0.0 # Because of vacuum BC
        for i,cell in enumerate(cells):
            psi_in = psi[i+1,n] = getPsiOut(psi_in,mu,cell,dx,method)

    # Backward Sweep
    for n in range(int(S.N/2),S.N):
        mu = S.mu[n]
        psi_in = 0.0 # Because of vacuum BC
        for i,cell in enumerate(cells):
            psi_in = psi[-2-i,n] = getPsiOut(psi_in,mu,cell,dx,method)


    for i,cell in enumerate(cells):
        flux_i = 0.0
        for n in range(S.N):
            wgtF = S.wgt[n]
            psiR = psi[i+1,n]
            psiL = psi[i,n]
            flux_i += wgtF*(psiR + psiL)*0.5

        cell.phi = flux_i



















