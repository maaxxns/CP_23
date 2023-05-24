import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt


data = np.genfromtxt('Psi.csv', unpack=True)
print(data.shape)
frame_number = 100

iterator = len(data[0,:]/frame_number)
x =  np.linspace(-10,10, int(2/0.01))

def update(n):
    plt.close()
    anidata = data[iterator*n, :]

    print(n/frame_number*100,"%")
    line.set_ydata(anidata)
    return line


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim(0,0.2)
ax.set_xlim(-10,10)
line, = ax.plot(x,data[:,0]) 
anim = FuncAnimation(fig,update,frames = frame_number, interval = 20)
anim.save('schrodinger.gif', fps=15)


