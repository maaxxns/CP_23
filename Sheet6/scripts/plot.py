import numpy as np 
import matplotlib.pyplot as plt 

data_r_28 = np.genfromtxt("bin/Lorentz_r_is_28.csv", comments="#", delimiter=", ") # data_r_28[:,0] this is X(t)
data_r_20 = np.genfromtxt("bin/Lorentz_r_is_20.csv", comments="#", delimiter=", ") # data_r_20[:,0] this is X(t)

########################################################

# projection

########################################################

plt.plot(data_r_20[:,0], data_r_20[:,1], lw=0.5)
plt.title("Projection onto xy plane")
plt.ylabel("y")
plt.xlabel("x") 
plt.savefig("build/Lorentz_r_20_projection.pdf")

plt.close()
plt.plot(data_r_28[:,0], data_r_28[:,1], lw=0.5)
plt.title("Projection onto xy plane")
plt.ylabel("y")
plt.xlabel("x")
plt.savefig("build/Lorentz_r_28_projection.pdf")

plt.close()

########################################################

# Poincare slice

########################################################

mask = np.ma.getmask(np.ma.masked_inside(data_r_20[:,2],19,20))

plt.scatter(data_r_20[:,0][mask], data_r_20[:,1][mask])
plt.title("Poincare slice of Z=20")
plt.ylabel("y")
plt.xlabel("x") 
plt.savefig("build/Lorentz_r_20_slice.pdf")

mask = np.ma.getmask(np.ma.masked_inside(data_r_28[:,2],19.5,20.5))

plt.close()
plt.scatter(data_r_28[:,0][mask], data_r_28[:,1][mask])
plt.title("Poincare slice of Z=20")
plt.ylabel("y")
plt.xlabel("x")
plt.savefig("build/Lorentz_r_28_slice.pdf")

plt.close()



#########################################################

# 3d plots 

#########################################################

ax = plt.figure().add_subplot(projection='3d')

ax.plot(*data_r_28.T, lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Lorenz Attractor r = 28")

plt.savefig("build/Lorentz_r_28_3d.pdf")

plt.close()



ax = plt.figure().add_subplot(projection='3d')

ax.plot(*data_r_20.T, lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Lorenz Attractor r = 20")

plt.savefig("build/Lorentz_r_20_3d.pdf")

plt.close()