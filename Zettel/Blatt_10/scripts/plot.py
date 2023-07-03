import numpy as np 
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt 

#get from: ./build/data.csv
#save in: ./build/Graph.pdf

H , m_num, m_ana = np.genfromtxt("./build/Ex1_data.csv", delimiter=' ', unpack= True)

plt.plot(H[::10],m_num[::10], '.',label='Numeric', alpha = 0.3)
plt.plot(H[::10],m_ana[::10], '.',label='Analytic', alpha = 0.3)
plt.legend()
plt.xlabel("H")
plt.ylabel("m")
plt.savefig("./build/Ex1_Graph.pdf")