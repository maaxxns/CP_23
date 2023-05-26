import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt


data = np.genfromtxt('Psi.csv')
#data = data.T
print(data.shape)
print(len(data[0,:]))
frame_number = 100

iterator = len(data[:,0]) // frame_number
x =  np.linspace(-10,10, data.shape[1])

def update(n):
    plt.close()
    anidata = data[iterator*n, :]

    print(f"\r{n/frame_number*100:.1f}%",end="")
    line.set_ydata(anidata)
    return line,


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim(0,4)
ax.set_xlim(-10,10)
line, = ax.plot(x,data[0,:]) 
anim = FuncAnimation(fig,update,frames = frame_number, interval = 20)
anim.save('schrodinger.gif', fps=15)
print()


