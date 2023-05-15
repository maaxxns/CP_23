import numpy as np



def Gauss_Seidel(Phi, Delta, epsilon, rho, start_condions):
    J = 1/Delta
    L = 1/Delta
    Phi_pre = Phi
    Phi[0,:] = start_condions[0]
    Phi[J,:] = start_condions[1]
    Phi[:,0] = start_condions[2]
    Phi[:,L] = start_condions[3]
    index = np.arange((1,J-1))
    while True:
        for j in index:
            for l in index:
                Phi[j,l] = 0.25 * (Phi[j+1,l] + Phi[j-1,l] + Phi[j,l+1] + Phi[j,l-1]) + 0.25 * Delta**2 * rho[j,l]
        if (np.norm(Phi-Phi_pre) < epsilon):
            break

    return Phi
    

