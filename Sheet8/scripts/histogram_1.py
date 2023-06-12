import numpy as np 
import matplotlib.pyplot as plt 

def histfunction(data, path):
    counts, bins = np.histogram(data, bins=10)
    path = "build/histogram_" + path
    plt.figure()
    plt.stairs(counts, bins)
    plt.xlabel("Value")
    plt.ylabel("counts")
    plt.title("histogram of " + path[16])
    plt.savefig(path)
    plt.close()

def plot_function(data, path):
    path = "build/plot_" + path
    plt.figure()
    plt.scatter(data[-1:-len(data)//2:-1], data[-2:len(data)//2-1:-1], marker=".")
    plt.xlabel(r"$r_n$")
    plt.ylabel(r"$r_{n + 1}$")
    plt.title("Coorelation of " + path[11])
    plt.savefig(path)
    plt.close()

 
a_data = np.genfromtxt("bin/rng_number_a.csv", comments="#")
b_data = np.genfromtxt("bin/rng_number_b.csv", comments="#")
c_data = np.genfromtxt("bin/rng_number_c.csv", comments="#")
d_data = np.genfromtxt("bin/rng_number_d.csv", comments="#")

datas = [a_data, b_data, c_data, d_data]
paths = ["a.pdf", "b.pdf", "c.pdf", "d.pdf"]
for i in range(len(datas)):
    histfunction(datas[i], paths[i])
    plot_function(datas[i], paths[i])