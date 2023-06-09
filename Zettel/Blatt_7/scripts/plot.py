import numpy as np 
import matplotlib.pyplot as plt 

#get from: ./build/data.csv
#save in: ./build/Graph.pdf


Ev_Householder = np.genfromtxt("./build/Ev_House_10.csv")
Ev_Potenz = np.genfromtxt("./build/Ev_Pot_10.csv")

Ev_Householder_sort = np.sort(Ev_Householder)
Ev_Potenz_sort = np.sort(Ev_Potenz)

x = np.arange(0,10)
plt.plot(x, Ev_Householder_sort, "x" ,label ='Householder')
plt.plot(x, Ev_Potenz_sort ,'x' ,label ='Potenzmethode')
plt.legend()
plt.savefig("./build/Ev_vergleich.pdf")
# Ex2

