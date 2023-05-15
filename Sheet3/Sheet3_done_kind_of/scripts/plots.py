import numpy as np 
import matplotlib.pyplot as plt 

N = str(256)

T = ["1", "1", "1"]
for i in range(len(T)):
    data = np.genfromtxt('bin/bN' + N + T[i] +'.dat', delimiter=',', comments="#")
    paar = np.genfromtxt('bin/bg_N' + N + T[i] + '.dat', delimiter=',', comments="#")

    plt.figure()
    plt.plot(data[1:,0], data[1:, 2], label='Ekin')
    plt.plot(data[1:,0], data[1:, 3], label='Epot') 
    plt.plot(data[1:,0], data[1:, 2] + data[1:, 3], label='Ekin + Epot')
    plt.xlabel("t")
    plt.ylabel("Energie")
    plt.legend()
    plt.title("Ekin, h = 0.01")
    plt.tight_layout()
    plt.savefig('plots/plot_E_' + N + T[i] + '.pdf')

    plt.close()

    plt.figure()
    plt.plot(data[:,0], data[:, 1], label='T')
    plt.xlabel("t")
    plt.ylabel("T")
    plt.legend()
    plt.title("Temperatur, h = 0.01")
    plt.tight_layout()
    plt.savefig('plots/plot_T_' + N + T[i] + '.pdf')

    plt.close()
    x = np.linspace(0,100,100)
    plt.figure()
    plt.bar(x,paar,label='Paircorreltion')
    plt.xlabel("r")
    plt.ylabel("counts")
    plt.legend()
    plt.title("Paircorreltion, h = 0.01")
    plt.tight_layout()
    plt.savefig('plots/paircorreclation' + N + T[i] + '.pdf')

    plt.close()
