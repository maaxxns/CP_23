import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

#get from: ./build/data.csv
#save in: ./build/Graph.pdf


data_global = np.genfromtxt("./build/data_global.csv", delimiter=' ', unpack=True)
data_perso = np.genfromtxt("./build/data_perso.csv", delimiter=' ', unpack=True)


iterations = np.arange(data_global.size)

for i in np.arange(data_perso.shape[0]):
    plt.plot(iterations, data_perso[i][:], '.', alpha = 0.4)
plt.plot(iterations, data_global, 'x', alpha = 1, label ='Global')
plt.legend()
plt.savefig("./build/Global_best.pdf")
plt.clf()


data_r = np.genfromtxt('./build/data_r.csv', delimiter= ' ', unpack=True)
data_v = np.genfromtxt('./build/data_v.csv', delimiter= ' ', unpack=True)
data_r_x = data_r[:][0:-2:2]
data_r_y = data_r[:][1:-1:2]
data_v_x = data_v[:][0:-2:2]
data_v_y = data_v[:][1:-1:2]
print(data_r_x.shape)
print(data_r_y.shape)
data_r = np.append(data_r_x,data_r_y).reshape(2,9,101).transpose()
print(data_r_x[:,80]) # first index are the x-cordiates and second the iteration

plt.scatter(data_r_x[:,0],data_r_y[:,0])
plt.savefig("./build/Test1.pdf")
plt.clf()

fig, ax = plt.subplots()
ax.set_xlim([-10, 10])
ax.set_ylim([-10, 10])

scat = ax.scatter(data_r_x, data_r_y)


def update(i):
    scat.set_offsets(data_r[:,:,i])
    return scat,

ani = FuncAnimation(fig, update, repeat=True, interval=10)

ani.save('./build/Test.gif', fps=15)




