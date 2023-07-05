import numpy as np 
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt 

#get from: ./build/data.csv
#save in: ./build/Graph.pdf

H , m_num, m_ana = np.genfromtxt("./build/Ex1_data.csv", delimiter=' ', unpack= True)

N = 30

plt.plot(H[::N],m_num[::N], '.',label='Numeric', alpha = 0.3)
plt.plot(H[::N],m_ana[::N], '.',label='Analytic', alpha = 0.3)
plt.legend()
plt.xlabel("H")
plt.ylabel("m")
plt.savefig("./build/Ex1_Graph.pdf")
plt.clf()

rel_abw = (m_num[::N]-m_ana[::N])/m_ana[::N]
H_abw = H[::N]
plt.plot(H_abw[rel_abw>-0.9],rel_abw[rel_abw>-0.9], '.',label='rel. Abweichung', alpha = 0.5)
plt.xlabel("H")
plt.ylabel("rel. Abweichung")
plt.savefig("./build/Ex1_Abweichung.pdf")
plt.clf()