import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

data = np.genfromtxt("bin/compressedwiggle_that_membrane.csv", comments="#", delimiter=", ") # data[0,:] gives array of time 0 so first index is time second is space kind of

iterator = len(data[:,0])//100

def animate(n):
    anidata = data[iterator*n,:].reshape(100,150)
    if(n == 0):
        print(0, "%")
    else:
        print(n, " %")
    extent = [0, 1.5, 0, 1]
    im = ax.imshow(anidata ,aspect="auto",vmin=-1, vmax=1, cmap="magma", extent=extent)
    return [im]

fig = plt.figure()
ax = fig.add_subplot(111)
anim = FuncAnimation(fig, animate, frames=100, interval=20, blit=True) #intervall is time between frames, blit is only frames that changed are redrawn
anim.save('build/animation.mp4', fps=15, extra_args=['-vcodec', 'libx264'])