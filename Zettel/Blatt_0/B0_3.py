import numpy as np
import matplotlib.pyplot as plt

#a)
N=100 #number of iterations



def euler(N,y_0,t_end):
    index = np.arange(N)
    del_t = t_end/N
    array = np.zeros(N+1)
    array[0] = y_0
    for x in index:
        array[x+1] = array[x]*(1-del_t)
    return array

def sym_euler(N,y_0, y_1,t_end):
    index = np.arange(1,N)
    del_t = t_end/N
    array = np.zeros(N+1)
    array[0] = y_0
    array[1] = y_1
    for x in index:
        array[x+1] = -2*del_t*array[x]+array[x-1]
    return array



t_end = 10
eu = euler(N,1,t_end)
sym_eu = sym_euler(N,1,np.exp(-t_end/N),t_end)
exakt = np.exp(-np.linspace(0,t_end,N+1))

t_lin = np.linspace(0,t_end,N+1)

plt.plot(t_lin,eu, label = 'euler')
plt.plot(t_lin,sym_eu, label = 'sym_euler')
plt.plot(t_lin,exakt, label = 'exakt')
plt.legend()
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("a.pdf")
plt.clf()


plt.plot(t_lin,eu, label = 'euler')
plt.plot(t_lin,sym_eu, label = 'sym_euler')
plt.plot(t_lin,exakt, label = 'exakt', ls= '--')
plt.legend()
plt.grid()
plt.xlim(0,2)
plt.ylim(0,1)
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("a_zoom.pdf")
plt.clf()
#b)

def euler_b(N,y_0,t_end):
    index = np.arange(N)
    del_t = t_end/N
    array = np.zeros(N+1)
    array[0] = y_0 - del_t
    for x in index:
        array[x+1] = array[x]*(1-del_t)
    return array


def sym_euler_b(N,y_0,t_end):
    index = np.arange(1,N)
    del_t = t_end/N
    array = np.zeros(N+1)
    array[0] = y_0
    array[1] = y_0-del_t
    for x in index:
        array[x+1] = -2*del_t*array[x]+array[x-1]
    return array



eu_b = euler_b(N,1,t_end)
sym_eu_b = sym_euler_b(N,1,t_end)


plt.plot(t_lin,eu_b, label = 'euler')
plt.plot(t_lin,sym_eu_b, label = 'sym_euler')
plt.plot(t_lin,exakt, label = 'exakt')
plt.legend()
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("b.pdf")
plt.clf()

plt.plot(t_lin,eu_b, label = 'euler')
plt.plot(t_lin,sym_eu_b, label = 'sym_euler')
plt.plot(t_lin,exakt, label = 'exakt', ls='--')
plt.legend()
plt.grid()
plt.xlim(0,2)
plt.ylim(0,1)
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("b_zoom.pdf")
plt.clf()