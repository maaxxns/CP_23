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



k_BT_array = np.array(["1.5", "2.25", "3"])
start_conditions_array = np.array(["ord" , "rnd"])
form_array = np.array(['x', '.', 'o'])
i = 0
for k_BT in k_BT_array:
    for start_condition in start_conditions_array:
        sweep, av_energy = np.genfromtxt("./build/av_energy_"+start_condition+"_k_BT_"+k_BT+".csv", unpack=True, delimiter=',')
        if sweep.size > 1000:
            plt.plot(sweep[::20], av_energy[::20], form_array[i], label = "k_BT="+k_BT+"_"+start_condition, alpha = 0.6)
        else:
            plt.plot(sweep, av_energy, form_array[i], label = "k_BT="+k_BT+"_"+start_condition, alpha = 0.6)
    i = i+1
plt.xlabel("sweep")
plt.ylabel("avg. Energy per Spin")
plt.title("Euqilibrate comparison")
plt.legend()
plt.tight_layout()
plt.savefig("./build/av_energy.pdf")
plt.clf()
