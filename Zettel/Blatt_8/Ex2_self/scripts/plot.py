import numpy as np 
import matplotlib.pyplot as plt 

#get from: ./build/data.csv 
#save in: ./build/Graph.pdf


Box = np.genfromtxt("./build/Box_data.csv", unpack=True)


Nr = np.array(["Box"])


for i in Nr:
    plt.hist(globals()[i],100, label = "Histogram " + i +")")
    plt.legend()
    plt.savefig("./build/" + i + ".pdf")
    plt.clf()


