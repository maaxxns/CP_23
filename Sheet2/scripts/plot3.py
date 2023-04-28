import numpy as np 
import matplotlib.pyplot as plt 

trajectory= np.genfromtxt('bin/euler.csv', delimiter=', ')
trajectory2= np.genfromtxt('bin/verlet.csv', delimiter=', ')
#print(trajectory)

plot_length = (len(trajectory)) # the .csv file contains the tracetory off normal euler/verlet + backwards euler/verlet. So we have to split the csv file in two part to get the results for back integration isolated

plt.figure()
plt.plot(trajectory[:plot_length//2,0], trajectory[:plot_length//2,1], label='mass1 euler') 
plt.plot(trajectory[:plot_length//2,2], trajectory[:plot_length//2,3], label='mass2 euler')
plt.plot(trajectory2[:plot_length//2,0], trajectory2[:plot_length//2,1], label='mass1 verlet') 
plt.plot(trajectory2[:plot_length//2,2], trajectory2[:plot_length//2,3], label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The trajectories of both masses, calculated with T = 100, h =0.001 \n with euler and verlet")
plt.tight_layout()
plt.savefig('plots/plota.pdf')

plt.close()


#print(trajectory)
plt.figure()
plt.plot(trajectory[:plot_length//2,0], trajectory[:plot_length//2,1],ls="-", label='mass1 euler') 
plt.plot(trajectory[:plot_length//2,2], trajectory[:plot_length//2,3],ls="-", label='mass2 euler')

plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The trajectories of both masses up to T = 100, h=0.001, \n calculated euler")
plt.tight_layout()
plt.savefig('plots/plotb.pdf')

plt.close()

plt.figure()
plt.plot(trajectory2[:plot_length//2,0], trajectory2[:plot_length//2,1],ls="-", label='mass1 verlet') 
plt.plot(trajectory2[:plot_length//2,2], trajectory2[:plot_length//2,3],ls="-", label='mass2 verlet')

plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The trajectories of both masses up to T = 100, h=0.001, \n calculated verlet")
plt.tight_layout()
plt.savefig('plots/plotc.pdf')

plt.close()

plt.figure()
plt.plot(trajectory[plot_length//2:,0], trajectory[plot_length//2:,1],ls="-", label='mass1 euler') 
plt.plot(trajectory[plot_length//2:,2], trajectory[plot_length//2:,3],ls="-", label='mass2 euler')
plt.plot(trajectory2[plot_length//2:,0], trajectory2[plot_length//2:,1],ls="-", label='mass1 verlet') 
plt.plot(trajectory2[plot_length//2:,2], trajectory2[plot_length//2:,3],ls="-", label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The back integrated trajectories of both masses with\n T = 100, h = -0.001, calculated with verlet and euler")
plt.tight_layout()
plt.savefig('plots/plotd.pdf')

plt.close()

plt.figure()
plt.plot(trajectory[plot_length//2:,0], trajectory[plot_length//2:,1],ls="-", label='mass1 euler') 
plt.plot(trajectory[plot_length//2:,2], trajectory[plot_length//2:,3],ls="-", label='mass2 euler')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The back integrated trajectories of both masses with\n T = 100, h = -0.001, calculated with euler")
plt.tight_layout()
plt.savefig('plots/plote.pdf')

plt.close()

plt.figure()
plt.plot(trajectory2[plot_length//2:,0], trajectory2[plot_length//2:,1],ls="-", label='mass1 verlet') 
plt.plot(trajectory2[plot_length//2:,2], trajectory2[plot_length//2:,3],ls="-", label='mass2 verlet')
plt.xlabel("x Position")
plt.ylabel("y Position")
plt.legend()
plt.title("The back integrated trajectories of both masses with\n T = 100, h = -0.001, calculated with verlet")
plt.tight_layout()
plt.savefig('plots/plotf.pdf')

plt.close()
