import numpy as np
import matplotlib.pyplot as plt

#Nr.1
print("Hello World")

#Nr.2
#a)

def a(x):
    return 1/np.sqrt(x)-1/(np.sqrt(x+1))


x = np.linspace(0,10**15,1000000)
print(a(x))
print(np.min(a(x)))
x_min = x[np.argmin(a(x))]
print(x_min)
print(a(x_min))

#b)
x = np.linspace(10**(-12),10,1000000)
def b(x):
    return (1-np.cos(x))/np.sin(x)



#c)
def c(x,delta):
    return np.sin(x+delta)-np.sin(x)

