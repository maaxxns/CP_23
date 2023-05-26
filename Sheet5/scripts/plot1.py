import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

frame_number = 1000


data = np.genfromtxt("bin/schroedinger.csv", comments="#", delimiter=", ") # data[0,:] gives array of time 0 so first index is time second is space kind of

#data[0,:] gives the vector at time 0
if(frame_number > len(data[:,0])):
    frame_number = len(data[:,0])
x = np.linspace(-10, 10, int(20/0.1))

iterator = len(data[:,0])//frame_number

def update(n):
    plt.close()
    anidata = data[iterator*n,:]
    if(n == 0):
        print(0, "%") 
    else:
        print(n/frame_number *100, " %")
    line.set_ydata(anidata)
    return line
    

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim(0, 0.2)
ax.set_xlim(-10,10)
ax.set_title('Probabilty distribution')
ax.set_xlabel("x")
ax.set_ylabel("Probability")
line, = ax.plot(x, data[0,:])
anim = FuncAnimation(fig, update, frames=frame_number, interval=20) #intervall is time between frames, blit is only frames that changed are redrawn
anim.save('build/schroedinger.gif',writer='pillow', fps=30) #, extra_args=['-vcodec', 'libx264']