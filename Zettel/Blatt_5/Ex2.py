import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import seaborn as sns


T_0 = np.genfromtxt("T_0.csv")
T_1 = np.genfromtxt("T_1.csv")

ax = sns.heatmap(T_0)
plt.show()
