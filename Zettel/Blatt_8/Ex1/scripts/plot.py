import numpy as np 
import matplotlib.pyplot as plt 

#get from: ./build/data.csv 
#save in: ./build/Graph.pdf


a = np.genfromtxt("./build/a_data.csv", unpack=True)
b = np.genfromtxt("./build/b_data.csv", unpack=True)
c = np.genfromtxt("./build/c_data.csv", unpack=True)
d = np.genfromtxt("./build/d_data.csv", unpack=True)


Nr = np.array(["a","b","c","d"])
Nr1 = np.array(["a","b"])
Nr2 = np.array(["c","d"])

for i in Nr:
    plt.hist(globals()[i],10, label = "Histogram " + i +")")
    plt.legend()
    plt.savefig("./build/" + i + ".pdf")
    plt.clf()



for i in Nr1:
    globals()[i+"_x"] = globals()[i][0:-2:2]
    globals()[i+"_y"] = globals()[i][1:-1:2]


    plt.scatter(globals()[i+"_x"],globals()[i+"_y"], alpha=0.4, s = 1)
    plt.title("Scatterplot of "+i+")")
    plt.savefig("./build/"+i+"2d.pdf")
    plt.clf()

for i in Nr2:
    globals()[i+"_x"] = globals()[i][0:-2:2]
    globals()[i+"_y"] = globals()[i][1:-1:2]


    plt.scatter(globals()[i+"_x"][::10],globals()[i+"_y"][::10], alpha=0.4, s = 1)
    plt.title("Scatterplot of "+i+")")
    plt.savefig("./build/"+i+"2d.pdf")
    plt.clf()

