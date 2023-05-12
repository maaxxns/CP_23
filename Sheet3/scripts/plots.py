import numpy as np 
import matplotlib.pyplot as plt 

data = np.genfromtxt('bin/b)set.dat', delimiter=',', comments="#")

plt.figure()
plt.plot(data[:,0], data[:, 2], label='Ekin')
plt.plot(data[:,0], data[:, 3], label='Epot') 
plt.plot(data[:,0], data[:, 2] + data[:, 3], label='Ekin + Epot')
plt.xlabel("t")
plt.ylabel("Energie")
plt.legend()
plt.title("Ekin, h = 0.01")
plt.tight_layout()
plt.savefig('plots/plot_E_16.pdf')

plt.close()

plt.figure()
plt.plot(data[:,0], data[:, 1], label='T')

plt.xlabel("t")
plt.ylabel("T")
plt.legend()
plt.title("Temperatur, h = 0.01")
plt.tight_layout()
plt.savefig('plots/plot_T_16.pdf')

plt.close()