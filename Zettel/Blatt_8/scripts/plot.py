import numpy as np 
import matplotlib.pyplot as plt 

#get from: ./build/data.csv 
#save in: ./build/Graph.pdf


a = np.genfromtxt("./build/a_data.csv", unpack=True)
b = np.genfromtxt("./build/b_data.csv", unpack=True)
c = np.genfromtxt("./build/c_data.csv", unpack=True)
d = np.genfromtxt("./build/d_data.csv", unpack=True)


Nr = np.array(["a","b","c","d"])

for i in Nr:
    plt.hist(globals()[i],10, label = "Histogram " + i +")")
    plt.legend()
    plt.savefig("./build/" + i + ".pdf")
    plt.clf()

