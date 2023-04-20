import numpy as np
import matplotlib.pyplot as plt

h, f = np.genfromtxt("bin/Ex_1_a.csv", delimiter=',', unpack=True)


plt.plot(h,f, 'x', label = 'numeric')
plt.axhline(y=np.cos(np.pi/4), color='r', linestyle='-',label = 'analytcal')
plt.xlabel('h')
plt.ylabel('f1´(x)')
plt.xscale('log')
plt.legend()
plt.grid()
plt.savefig('bin/Ex_1_a.pdf')
plt.clf()


x1, del_f1 = np.genfromtxt("bin/Ex_1_a_2.csv", delimiter=',', unpack=True)

cut = (-0.05<del_f1) & (del_f1<0.05)
plt.plot(x1,del_f1, '.', label = 'rel_error')
plt.xlabel('x')
plt.ylabel('f1(x)')
plt.legend()
plt.grid()
plt.savefig('bin/Ex_1_a_2.pdf')
plt.clf()
plt.plot(x1[cut],del_f1[cut], '.', label = 'rel_error_cut')
plt.xlabel('x')
plt.ylabel('f1(x)')
plt.legend()
plt.grid()
plt.savefig('bin/Ex_1_a_2_cut.pdf')
plt.clf()


######b)

h, f = np.genfromtxt("bin/Ex_1_b.csv", delimiter=',', unpack=True)


plt.plot(h,f, 'x', label = 'numeric')
plt.axhline(y=np.cos(np.pi/4), color='r', linestyle='-',label = 'analytcal')
plt.xlabel('h')
plt.ylabel('f1´(x)')
plt.xscale('log')
plt.legend()
plt.grid()
plt.savefig('bin/Ex_1_b.pdf')
plt.clf()

x, del_f = np.genfromtxt("bin/Ex_1_b_2.csv", delimiter=',', unpack=True)

plt.plot(x,del_f, '.', label = 'rel_error')
plt.xlabel('x')
plt.ylabel('f1(x)')
plt.legend()
plt.grid()
plt.savefig('bin/Ex_1_b_2.pdf')
plt.clf()


#####c)
x, del_f = np.genfromtxt("bin/Ex_1_c.csv", delimiter=',', unpack=True)
cut = (-0.05<del_f1) & (del_f1<0.05)


plt.plot(x1[cut],del_f1[cut], '.', label = 'rel_error_two')
plt.plot(x[cut],del_f[cut], '.', label = 'rel_error_four')
plt.xlabel('x')
plt.ylabel('f1(x)')
plt.legend()
plt.grid()
plt.savefig('bin/Ex_1_c.pdf')
plt.clf()