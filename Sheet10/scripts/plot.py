import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

data = np.genfromtxt("bin/state.csv", comments="#", delimiter=", ") # data[0,:] gives array of time 0 so first index is time second is space kind of

data_1 = np.genfromtxt("bin/state_at_T_1.csv", delimiter=", ")
data_3 = np.genfromtxt("bin/state_at_T_3.csv", delimiter=", ")

data_energy_1_5 = np.genfromtxt("bin/energy_per_spinT_1_5.csv", delimiter=", ")
data_energy_3 = np.genfromtxt("bin/energy_per_spinT_3.csv", delimiter=", ")

data_abs_mag_1_5 = np.genfromtxt("bin/abs_magnetT_1_5.csv", delimiter=", ")
data_abs_mag_3 = np.genfromtxt("bin/abs_magnetT_3.csv", delimiter=", ")

data_mag_1_5 = np.genfromtxt("bin/magnetT_1_5.csv", delimiter=", ")
data_mag_3 = np.genfromtxt("bin/magnetT_3.csv", delimiter=", ")


extent = [0, 100, 0, 100]
plt.figure()
plt.imshow(data_1.reshape(100, 100), aspect="auto",vmin=-1, vmax=1, cmap="magma", extent=extent)
plt.title("Snapshot at kT = 1")
plt.savefig("build/snapshotT_1.pdf")
plt.close()

plt.figure()
extent = [0, 100, 0, 100]
plt.imshow(data_3.reshape(100, 100), aspect="auto",vmin=-1, vmax=1, cmap="magma", extent=extent)
plt.title("Snapshot at kT = 3")
plt.savefig("build/snapshotT_3.pdf")
plt.close()
###################################################################################################
x = np.linspace(0, len(data_energy_1_5), len(data_energy_1_5))

plt.figure()
plt.plot(x, data_energy_1_5, "ko", label="Energy per State")
plt.title("Energy per State per Time step for 1.5=Tk")
plt.xlabel("time/steps")
plt.ylabel("Energy")
plt.savefig("build/energy_per_state_T_1_5.pdf")
plt.close()

plt.figure()
plt.plot(x, data_energy_3, "ko", label="Energy per State")
plt.title("Energy per State per Time step for 3=Tk")
plt.xlabel("time/steps")
plt.ylabel("Energy")
plt.savefig("build/energy_per_state_T_3.pdf")
plt.close()

###################################################################################################
plt.figure()
plt.plot(x, data_mag_1_5, "ko", label="magnetization")
plt.title("Magnetization for 1.5=Tk")
plt.xlabel("time/steps")
plt.ylabel(r"$\sum m$")
plt.savefig("build/magnet_T_1_5.pdf")
plt.close()

plt.figure()
plt.plot(x, data_mag_3, "ko", label="magnetization")
plt.title("Magnetization for 3=Tk")
plt.xlabel("time/steps")
plt.ylabel(r"$\sum m$")
plt.savefig("build/magnet_T_3.pdf")
plt.close()

###################################################################################################

plt.figure()
plt.plot(x, data_abs_mag_1_5, "ko", label="magnetization")
plt.title("Magnetization for 1.5=Tk")
plt.xlabel("time/steps")
plt.ylabel(r"$\sum |m|$")
plt.savefig("build/abs_magnet_T_1_5.pdf")
plt.close()

plt.figure()
plt.plot(x, data_abs_mag_3, "ko", label="magnetization")
plt.title("Magnetization for 3=Tk")
plt.xlabel("time/steps")
plt.ylabel(r"$\sum |m|$")
plt.savefig("build/abs_magnet_T_3.pdf")
plt.close()



###################################################################################################

iterator = len(data[:,0])//100

def animate(n):
    anidata = data[iterator*n,:].reshape(100,100)
    if(n == 0):
        print(0, "%\r")
    else:
        print(n, " %\r")
    extent = [0, 100, 0, 100]
    im = ax.imshow(anidata ,aspect="auto",vmin=-1, vmax=1, cmap="magma", extent=extent)
    return [im]



fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Ising')
ax.set_xlabel("x")
ax.set_ylabel("y")
anim = FuncAnimation(fig, animate, frames=100, interval=20, blit=True) #intervall is time between frames, blit is only frames that changed are redrawn
anim.save('build/ising.gif',writer='pillow', fps=15) #, extra_args=['-vcodec', 'libx264']