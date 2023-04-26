import numpy as np 
import matplotlib.pyplot as plt 

trajectory= np.genfromtxt('src/euler.csv', delimiter=', ')
trajectory2= np.genfromtxt('src/verlet.csv', delimiter=', ')
#print(trajectory)

plt.figure()
plt.plot(trajectory[:1000,0], trajectory[:1000,1], label='mass1 euler') 
plt.plot(trajectory[:1000,2], trajectory[:1000,3], label='mass2 euler')
plt.plot(trajectory2[:1000,0], trajectory2[:1000,1], label='mass1 verlet') 
plt.plot(trajectory2[:1000,2], trajectory2[:1000,3], label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The trajectories of both masses, calculated with euler and verlet")
plt.tight_layout()
plt.savefig('plots/plota.pdf')

plt.close()


#print(trajectory)
plt.figure()
plt.plot(trajectory[:100,0], trajectory[:100,1],ls="-", label='mass1 euler') 
plt.plot(trajectory[:100,2], trajectory[:100,3],ls="-", label='mass2 euler')
plt.plot(trajectory2[:100,0], trajectory2[:100,1],ls="-", label='mass1 verlet') 
plt.plot(trajectory2[:100,2], trajectory2[:100,3],ls="-", label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The trajectories of both masses up to T = 10 , calculated with verlet and euler")
plt.tight_layout()
plt.savefig('plots/plotb.pdf')

plt.close()

plt.figure()
plt.plot(trajectory[:,0], trajectory[:,1],ls="-", label='mass1 euler') 
plt.plot(trajectory[:,2], trajectory[:,3],ls="-", label='mass2 euler')
plt.plot(trajectory2[:,0], trajectory2[:,1],ls="-", label='mass1 verlet') 
plt.plot(trajectory2[:,2], trajectory2[:,3],ls="-", label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The trajectories of both masses with back integration -h, calculated with verlet and euler")
plt.tight_layout()
plt.savefig('plots/plotc.pdf')

plt.close()

plt.figure()
plt.plot(trajectory[1000:,0], trajectory[1000:,1],ls="-", label='mass1 euler') 
plt.plot(trajectory[1000:,2], trajectory[1000:,3],ls="-", label='mass2 euler')
plt.plot(trajectory2[1000:,0], trajectory2[1000:,1],ls="-", label='mass1 verlet') 
plt.plot(trajectory2[1000:,2], trajectory2[1000:,3],ls="-", label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("Just the back integration -h, calculated with verlet and euler")
plt.tight_layout()
plt.savefig('plots/plotd.pdf')

plt.close()

plt.figure()
plt.plot(trajectory[1000:1100,0], trajectory[1000:1100,1],ls="-", label='mass1 euler') 
plt.plot(trajectory[1000:1100,2], trajectory[1000:1100,3],ls="-", label='mass2 euler')
plt.plot(trajectory2[1000:1100,0], trajectory2[1000:1100,1],ls="-", label='mass1 verlet') 
plt.plot(trajectory2[1000:1100,2], trajectory2[1000:1100,3],ls="-", label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("Just the back integration from T=100-90 -h, calculated with verlet and euler")
plt.tight_layout()
plt.savefig('plots/plote.pdf')
