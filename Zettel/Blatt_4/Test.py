import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
#Delta = 0.5
#J = int(1/Delta)
#L = int(1/Delta)
#x_array = np.arange(0,J,Delta)
#y_array = np.arange(0,L,Delta)
#
#print(x_array)
#x, y = np.meshgrid(x_array, y_array)
#
#result = x+y
#print(x)
#print(result)





def gradient_operator(potential):
    # Calculate the gradient using central differences
    gradient_x = np.gradient(potential, axis=0)
    gradient_y = np.gradient(potential, axis=1)

    # Return the gradient components
    return gradient_x, gradient_y

# Example 2D potential
potential = np.array([[1, 2, 3],
                      [4, 5, 6],
                      [7, 8, 9]])

# Calculate the gradient
gradient_x, gradient_y = gradient_operator(potential)

# Print the gradient components
print("Gradient X:", gradient_x)
print("Gradient Y:", gradient_y)

plt.quiver([0,1,2], [0,1,2], gradient_x, gradient_y)
plt.show()
