




def getQ2(cell):
    return cell.src + 0.5*cell.SigS*cell.phi


def diamondDiff2(psi_in,mu,cell,dx):
    Q = getQ2(cell)
    return ( psi_in * (2*mu - cell.SigT*dx) + 2.0*Q*dx) / \
              ( 2*mu + cell.SigT*dx )






def transportSweep(quad,cells,dx):

    # Forward Sweep
    mu = abs(quad.mu[0])
    psi_forward = [0.0]*(len(cells)+1)
    psi_in = 0.0
    for i,cell in enumerate(cells):
        #Q = 0.5*(cell.SigS*cell.phi + cell.src)
        #psi_out = ( psi_in * (2.0*mu - dx*cell.SigT) + 2.0*dx*Q ) / \
        #          ( 2.0*mu + dx*cell.SigT )
        psi_out = diamondDiff2(psi_in,mu,cell,dx)
        #print(psi_out)
        psi_forward[i+1] = psi_out
        psi_in = psi_out
        #print("-----",psi_in)
    print(psi_forward)


    # Backward Sweep
    mu = abs(quad.mu[1])
    psi_backward = [0.0]*(len(cells)+1)
    psi_in = 0.0
    for i,cell in enumerate(cells):
        Q = 0.5*(cell.SigS*cell.phi + cell.src)
        psi_out = (psi_in*(2.0*mu - dx*cell.SigT) + 2.0*dx*Q) / \
                  (2.0 *mu + dx*cell.SigT)
        psi_backward[-2-i] = psi_out
        psi_in = psi_out

    psi_backward[0] = psi_out

    print(psi_backward)

    for i,cell in enumerate(cells):

        wgtF = quad.wgt[0]
        psiR = psi_forward[i+1]
        psiL = psi_forward[i]
        flux_i = wgtF*(psiR + psiL)/2.0

        wgtF = quad.wgt[1]
        psiR = psi_backward[i+1]
        psiL = psi_backward[i]
        flux_i += wgtF*(psiR + psiL)/2.0

        cell.phi = flux_i




















