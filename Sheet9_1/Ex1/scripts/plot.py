import numpy as np 
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt 

#get from: ./build/data.csv
#save in: ./build/Graph.pdf


iterations_h, x , y, z = np.genfromtxt("./build/halbierung.csv", delimiter=' ', unpack=True)

iterations_N, x_N = np.genfromtxt("./build/newton.csv", delimiter=' ', unpack=True)




plt.plot(iterations_h,x,'x' ,label ='x_i')
plt.plot(iterations_h,y,'x',label ='y_i')
plt.plot(iterations_h,z,'x',label ='z_i')
plt.xlabel("Iteration")
plt.legend()
plt.savefig("./build/Halbierung.pdf")
plt.clf()

plt.plot(iterations_N,x_N,'x' ,label ='x, Newton')
plt.xlabel("Iteration")
plt.legend()
plt.savefig("./build/Newton.pdf")
plt.clf()