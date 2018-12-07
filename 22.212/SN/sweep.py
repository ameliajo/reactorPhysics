import numpy as np




def getQ(cell):
    return cell.src + 0.5*cell.SigS*cell.phi


def step(psi_in,mu_n,cell):
    Q = getQ(cell)
    mu = abs(mu_n)
    return (cell.dx*Q/mu) + psi_in * (1.0 - cell.dx*cell.SigT/mu)


def diamondDiff(psi_in,mu_n,cell):
    Q = getQ(cell)
    mu = abs(mu_n)
    return ( ( 2.0*mu - cell.SigT*cell.dx ) * psi_in + 2.0*Q*cell.dx ) / \
             ( 2.0*mu + cell.SigT*cell.dx )

def stepCharacteristic(psi_in,mu_n,cell):
    Q = getQ(cell)
    mu = abs(mu_n)
    return (1.0 - (cell.dx*cell.SigT/mu)*np.exp(-cell.SigT*cell.dx/mu)) * psi_in + \
           ((cell.dx*Q/mu)*np.exp(-cell.SigT*cell.dx/mu))


       
def linearDiscontinuous(psi_in,mu_n,cell):
    Q = getQ(cell)
    #mu = abs(mu_n)
    mu = mu_n

    Q_R = Q
    Q_L = Q

    Q_tilde_L = (1.0*Q_R/3.0 + 2.0*Q_L/3.0)*cell.dx
    Q_tilde_R = (2.0*Q_R/3.0 + 1.0*Q_L/3.0)*cell.dx

    if mu > 0.0:
        T_23 = mu + 2.0*cell.SigT*cell.dx/3.0
        S_13 = mu - 1.0*cell.SigT*cell.dx/3.0

        # I'm given psi_R_i-1 and I'm trying to get psi_L_i
        psi_L_i = (Q_tilde_L + 2.0*mu*psi_in - cell.SigT*cell.dx*Q_tilde_R/(3.0*T_23) - mu*Q_tilde_R/T_23) /\
                  (mu + mu*S_13/T_23 + 2.0*cell.SigT*cell.dx/3.0 + cell.SigT*cell.dx*S_13/(3.0*T_23))
        
        # Now I have psi_L_i and I'm trying to get psi_R_i
        return Q_tilde_R/T_23 + psi_L_i*S_13/T_23

    if mu < 0.0:
        T_13 = mu + 1.0*cell.SigT*cell.dx/3.0
        S_23 = mu - 2.0*cell.SigT*cell.dx/3.0

        # I'm given psi_L_i+1 and I'm trying to get psi_R_i
        psi_R_i = (Q_tilde_R - 2.0*mu*psi_in + cell.SigT*cell.dx*Q_tilde_L/(3.0*S_23) - mu*Q_tilde_L/S_23) /\
                  (-mu -mu*T_13/S_23 + 2.0*cell.SigT*cell.dx/3.0 + cell.SigT*cell.dx*T_13/(3.0*S_23))

        # Now I have psi_R_i and I'm trying to get psi_L_i
        return psi_R_i*T_13/S_23 - Q_tilde_L/S_23






def getPsiOut(psi_in,mu,cell,method):
    if method == 'step':                return step(psi_in,mu,cell)
    if method == 'diamond':             return diamondDiff(psi_in,mu,cell)
    if method == 'stepCharacteristic':  return stepCharacteristic(psi_in,mu,cell)
    if method == 'linearDiscontinuous': return linearDiscontinuous(psi_in,mu,cell)
    raise ValueError

 



def broom(S,cells,method):
    assert(S.N%2 == 0) # Require that order is even

    psi = np.zeros((len(cells)+1,S.N))

    # Forward Sweep
    for n in range(int(S.N/2)):
        mu = S.mu[n]
        psi_in = 0.0 # Because of vacuum BC
        for i,cell in enumerate(cells):
            psi_in = psi[i+1,n] = getPsiOut(psi_in,mu,cell,method)

    # Backward Sweep
    for n in range(int(S.N/2),S.N):
        mu = S.mu[n]
        psi_in = 0.0 # Because of vacuum BC
        for i,cell in enumerate(cells):
            psi_in = psi[-2-i,n] = getPsiOut(psi_in,mu,cell,method)


    for i,cell in enumerate(cells):
        flux_i = 0.0
        for n in range(S.N):
            wgtF = S.wgt[n]
            psiR = psi[i+1,n]
            psiL = psi[i,n]
            flux_i += wgtF*(psiR + psiL)*0.5

        cell.phi = flux_i



















