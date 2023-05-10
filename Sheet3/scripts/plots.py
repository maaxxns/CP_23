import numpy as np 
import matplotlib.pyplot as plt 

data = np.genfromtxt('bin/b)set.dat', delimiter=',', comments="#")

plt.figure()
plt.plot(data[:,0], data[:, 2], label='Ekin') 
plt.xlabel("t")
plt.ylabel("Ekin")
plt.legend()
plt.title("Ekin, h = 0.01")
plt.tight_layout()
plt.savefig('plots/plot.pdf')

plt.close()