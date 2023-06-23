import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

test_x = np.array([[1,1],[2,1],[3,2],[4,3],[5,4],[6,5]])
test_y = np.array([[8,2],[9,3],[10,4],[11,5],[12,6],[13,7]])
test = np.append(test_x,test_y).reshape(2,6,2).transpose()
print(test)
