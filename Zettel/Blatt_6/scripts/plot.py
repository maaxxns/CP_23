import numpy as np 
import matplotlib.pyplot as plt 


x_0 = np.array([0.1,0.3,0.5,0.7,0.9])
for i in np.arange(0,len(x_0)):
    r , x = np.genfromtxt("./build/r_data_log"+str(i)+".csv", delimiter=',', unpack=True)

    
    

    plt.plot(r,x, '.', label = "x_0="+str(x_0[i]))
    plt.xlabel('r')
    plt.ylabel('Fixpunkt')
    plt.legend()
    plt.grid()
    plt.savefig("./build/Bif_log"+str(i)+".pdf")
    plt.clf()

for i in np.arange(0,len(x_0)):
    r , x = np.genfromtxt("./build/r_data_kub"+str(i)+".csv", delimiter=',', unpack=True)

    
    

    plt.plot(r,x, '.', label = "x_0="+str(x_0[i]))
    plt.xlabel('r')
    plt.ylabel('Fixpunkt')
    plt.legend()
    plt.grid()
    plt.savefig("./build/Bif_kub"+str(i)+".pdf")
    plt.clf()