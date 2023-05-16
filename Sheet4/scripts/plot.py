import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

data1 = np.genfromtxt("bin/Diffusion_a.csv", comments="#", delimiter=", ")
data2 = np.genfromtxt("bin/Diffusion_b_bad.csv", comments="#", delimiter=", ")
data3 = np.genfromtxt("bin/Diffusion_b_good.csv", comments="#", delimiter=", ")
data4 = np.genfromtxt("bin/Diffusion_c_haevy.csv", comments="#", delimiter=", ")
data5 = np.genfromtxt("bin/Diffusion_c_weird.csv", comments="#", delimiter=", ")

data = [data1, data2, data3, data4, data5]
filenames = ["2a.mp4", "2b_bad.mp4", "2b_good.mp4", "2c_haevy.mp4", "2c_weird.mp4"]

weird = False

for i in range(len(data)):
    print("processing file :", filenames[i])
    iterator = len(data[i])//100
    def animate(n):
        if(n == 0):
            print(0, "%")
        else:
            print(n, " %")
        #print(data[i][:,i*100])
        hist = data[i][:,n*iterator][np.newaxis, :]
        extent = [0, 100, 0, 1]
        if(i == 3):
            im = ax.imshow(hist ,aspect="auto", vmin=0, vmax=0.02, cmap="magma", extent=extent)
            return [im]
        if(i == 4):
            im = ax.imshow(hist ,aspect="auto", vmin=0, vmax=0.11, cmap="magma", extent=extent)
            return [im]
        else:
            im = ax.imshow(hist ,aspect="auto", vmin=0, vmax=1, cmap="magma", extent=extent)
            return [im]


    fig = plt.figure()
    ax = fig.add_subplot(111)

    anim = FuncAnimation(fig, animate, frames=100, interval=20, blit=True) #intervall is time between frames, blit is only frames that changed are redrawn

    anim.save('build/' + filenames[i], fps=15, extra_args=['-vcodec', 'libx264'])


