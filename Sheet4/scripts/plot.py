import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

data = np.genfromtxt("bin/Diffusion_a.csv", comments="#", delimiter=", ")

weird = True

print(len(data[0,:]))

def animate(i):
    if(i == 0):
        print(0)
    else:
        print(i*1000/len(data[0,:]) * 100)
    #print(data[:,i*100])
    hist = data[:,i*1000][np.newaxis, :]
    extent = [0, 100, 0, 1]
    if(weird == True):
        im = ax.imshow(hist ,aspect="auto", vmin=0, vmax=0.11, cmap="magma", extent=extent)
    else:
        im = ax.imshow(hist ,aspect="auto", vmin=0, vmax=1, cmap="magma", extent=extent)
    return [im]


fig = plt.figure()
ax = fig.add_subplot(111)

anim = FuncAnimation(fig, animate, frames=len(data[0,:])//1000, interval=20, blit=True) #intervall is time between frames, blit is only frames that changed are redrawn

anim.save('build/2a.mp4', fps=15, extra_args=['-vcodec', 'libx264'])


