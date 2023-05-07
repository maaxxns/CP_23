import numpy as np
import matplotlib.pyplot as plt


r1_x, r1_y, v1_x, v1_y = np.genfromtxt('bin/Ex_3_euler_1.csv' ,delimiter=',', unpack=True)
r2_x, r2_y, v2_x, v2_y = np.genfromtxt('bin/Ex_3_euler_2.csv' ,delimiter=',', unpack=True)

plt.plot(r1_x, r1_y, 'x')
plt.plot(r2_x, r2_y, 'x')
plt.savefig('build/test.pdf')
