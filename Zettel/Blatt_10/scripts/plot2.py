import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

#get from: ./build/data.csv
#save in: ./build/Graph.pdf


def create_color_map_from_csv(k_BT, Nr):
    data = np.loadtxt("./build/k_BT:"+str(k_BT)+"_sweep:"+str(Nr)+".csv", delimiter=',')
    cmap = mcolors.ListedColormap(['blue', 'red'])
    plt.imshow(data, cmap=cmap)
    plt.colorbar(ticks=[-1,1])
    plt.savefig("./build/Ising"+str(k_BT)+"_"+str(Nr)+".pdf")
    plt.title("k_BT:"+str(k_BT)+", Sweep:"+str(Nr*2e3))
    plt.tight_layout()
    plt.clf()

for k in np.array([1,3]):
    for i in np.array([0,1,2,3,4]):
        create_color_map_from_csv(k,i)



k_BT_array = np.array([1.5, 2.0, 2.25, 2.5, 3.0])
start_conditions_array = np.array("orderly" , "rnd")

for start_condiotion in start_conditions_array:
    for k_BT in k_BT_array:
        sweep, av_energy = np.genfromtxt("./build/av_energy_"+start_condiotion+"_k_BT_"+k_BT+".csv")
        plt.plot(sweep, av_energy, label = )
