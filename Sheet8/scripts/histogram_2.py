import numpy as np 
import matplotlib.pyplot as plt 

def histfunction(data, path):
    counts, bins = np.histogram(data, bins=100)
    plt.figure()
    plt.stairs(counts, bins)
    plt.xlabel("Value")
    plt.ylabel("counts")
    plt.title("histogram of muller")
    plt.savefig(path)
    plt.close()

def p(x):
    return np.sin(x)/2

def p_2(x):
    return 3 * x**2

data_muller = np.genfromtxt("bin/2_box_muller.csv", comments="#")
data_rejection = np.genfromtxt("bin/2_rejection.csv", delimiter=", ")
data_inverse = np.genfromtxt("bin/2_inverse.csv", delimiter=", ")



histfunction(data_muller, "build/muller_hist.pdf")

x = np.linspace(0, np.pi, 1000)

plt.figure()
plt.scatter(data_rejection[:,0], data_rejection[:,1],marker=".", label="drawn numbers")
plt.plot(x, p(x),c="k", ls="-" , label="p(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.title("The rejection scheme with p(x)")
plt.savefig("build/rejection.pdf")


plt.figure()
counts, bins = np.histogram(data_inverse, bins=100)
plt.stairs(counts, bins)
plt.plot(x/np.pi, p_2(x/np.pi),c="k", ls="-" , label="p2(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.title("The inverse scheme with p_2(x)")
plt.savefig("build/inverse.pdf")
