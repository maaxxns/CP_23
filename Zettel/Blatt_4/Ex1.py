import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def n_grad(potential):
    # Calculate the gradient using central differences
    gradient_x = np.gradient(potential, axis=1)
    gradient_y = np.gradient(potential, axis=0)

    # Return the negative gradient components
    return -gradient_x, -gradient_y


def Gauss_Seidel(Phi, Delta, epsilon, rho, start_condions):
    J = int(1/Delta)
    L = int(1/Delta)
    Phi[0,:] = start_condions[0]
    Phi[J,:] = start_condions[1]
    Phi[:,0] = start_condions[2]
    Phi[:,L] = start_condions[3]
    
    index = np.arange(1,J)
    count = 0
    while True:
        Phi_pre = Phi.copy()
        count += 1
        for j in index:
            for l in index:
                Phi[j,l] = 0.25 * (Phi[j+1,l] + Phi[j-1,l] + Phi[j,l+1] + Phi[j,l-1]) + 0.25 * Delta**2 * rho[j,l]
        if ((np.abs(Phi-Phi_pre) < epsilon).all()):
            break
    print("Iterations:", count)
    E_x , E_y = n_grad(Phi)
    return Phi , E_x , E_y
    

Delta = 0.05
epsilon = 1e-5


#b)
J = int(1/Delta)
L = int(1/Delta)

Phi_b = np.ones((J+1,J+1))
rho_b = np.zeros((J+1,J+1))
start_conditions_b = np.zeros(4)

Phi_solution_b = Gauss_Seidel(Phi_b, Delta, epsilon, rho_b, start_conditions_b)
ax = sns.heatmap(Phi_solution_b[0])
plt.savefig("build/Phi_b.pdf")
plt.clf()



#c)
Phi_c = np.ones((J+1,J+1))
rho_c = np.zeros((J+1,J+1))
start_conditions_c = np.array([0,1,0,0])
Phi_solution_c = Gauss_Seidel(Phi_c, Delta, epsilon, rho_c, start_conditions_c)
ax = sns.heatmap(Phi_solution_c[0])
plt.savefig("build/Phi_c.pdf")
plt.clf()

### analytic solution:
def analytic_solution(x_values, y_values, n_iterations):
    x, y = np.meshgrid(x_values, y_values)
    n_index = np.arange(1, n_iterations)
    result = np.zeros_like(x)
    for n in n_index:
        result += 2 * (1-np.cos(n*np.pi)) / (n*np.pi*np.sinh(n*np.pi)) * np.sin(n*np.pi*x) * np.sinh(n*np.pi*y)
    return result

x_array = np.arange(0,1+Delta,Delta)
y_array = np.arange(0,1+Delta,Delta)
Phi_analytic = analytic_solution(x_array,y_array,200)



ax = sns.heatmap(Phi_analytic)
plt.savefig("build/Phi_analytic.pdf")
plt.clf()
ax = sns.heatmap(np.abs(Phi_analytic-Phi_solution_c[0]))
plt.savefig("build/Phi_difference.pdf")
plt.clf()


#d)
Phi_d = np.zeros((J+1,J+1))
rho_d = np.zeros((J+1,J+1))
rho_d[int(J/2),int(L/2)] = 1
print(rho_d)
start_conditions_d = np.zeros(4)
Phi_solution_d = Gauss_Seidel(Phi_d, Delta, epsilon, rho_d, start_conditions_d)
ax = sns.heatmap(Phi_solution_d[0])
plt.savefig("build/Phi_d.pdf")
plt.clf()


plt.quiver(x_array, y_array, Phi_solution_d[1], Phi_solution_d[2] )
plt.savefig("build/E_d.pdf")
plt.clf()


#e)

Phi_e = np.zeros((J+1,J+1))
rho_e = np.zeros((J+1,J+1))
rho_e[int(J/4),int(L/4)] = 1
rho_e[int(3*J/4),int(3*L/4)] = 1
rho_e[int(J/4),int(3*L/4)] = -1
rho_e[int(3*J/4),int(L/4)] = -1
start_conditions_e = np.zeros(4)
Phi_solution_e = Gauss_Seidel(Phi_e, Delta, epsilon, rho_e, start_conditions_e)
ax = sns.heatmap(Phi_solution_e[0])
plt.savefig("build/Phi_e.pdf")
plt.clf()

plt.quiver(x_array, y_array, Phi_solution_e[1], Phi_solution_e[2] )
plt.savefig("build/E_e.pdf")
plt.clf()

