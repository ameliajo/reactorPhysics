




def getQ(cell):
    return cell.src + 0.5*cell.SigS*cell.phi


def diamondDiff(psi_in,mu,cell,dx):
    Q = getQ(cell)
    return ( psi_in * (2*mu - cell.SigT*dx) + 2.0*Q*dx) / \
              ( 2*mu + cell.SigT*dx )






def broom(quad,cells,dx):
    # Forward Sweep
    mu = abs(quad.mu[0])
    psi_forward = [0.0]*(len(cells)+1)
    psi_in = 0.0 # Because of vacuum BC
    for i,cell in enumerate(cells):
        psi_out = diamondDiff(psi_in,mu,cell,dx)
        psi_forward[i+1] = psi_out
        psi_in = psi_out

    # Backward Sweep
    mu = abs(quad.mu[1])
    psi_backward = [0.0]*(len(cells)+1)
    psi_in = 0.0 # Because of vacuum BC
    for i,cell in enumerate(cells):
        psi_out = diamondDiff(psi_in,mu,cell,dx)
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




















