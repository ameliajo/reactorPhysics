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
    if method == 'step':
        return step(psi_in,mu,cell,dx)
    if method == 'diamond':
        return diamondDiff(psi_in,mu,cell,dx)
    if method == 'step-characteristic':
        return stepCharacteristic(psi_in,mu,cell,dx)
    raise ValueError

        


def broom(quad,cells,dx,method):
    # Forward Sweep
    mu = quad.mu[0]
    psi_forward = [0.0]*(len(cells)+1)
    psi_in = 0.0 # Because of vacuum BC
    for i,cell in enumerate(cells):
        psi_out = getPsiOut(psi_in,mu,cell,dx,method)
        psi_forward[i+1] = psi_out
        psi_in = psi_out

    # Backward Sweep
    mu = quad.mu[1]
    psi_backward = [0.0]*(len(cells)+1)
    psi_in = 0.0 # Because of vacuum BC
    for i,cell in enumerate(cells):
        psi_out = getPsiOut(psi_in,mu,cell,dx,method)
        psi_backward[-2-i] = psi_out
        psi_in = psi_out


    for i,cell in enumerate(cells):

        wgtF = quad.wgt[0]
        psiR = psi_forward[i+1]
        psiL = psi_forward[i]
        flux_i = wgtF*(psiR + psiL)*0.5

        wgtF = quad.wgt[1]
        psiR = psi_backward[i+1]
        psiL = psi_backward[i]
        flux_i += wgtF*(psiR + psiL)*0.5

        cell.phi = flux_i




















