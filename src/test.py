import numpy as np
import matplotlib.pyplot as plt


tau = 12
#space = np.linspace(1/(tau+1), 1, tau+1)
#space = np.logspace(1/(tau+1), 1, tau+1, base=10)
space = np.full(tau+1,1)
print(space)
# Calculate the sum of the linspace
sum_linspace = np.sum(space)
print(sum_linspace)
# Normalize the linspace to make its sum equal to 1
normalized_linspace = space / sum_linspace
print(normalized_linspace)
print(np.sum(normalized_linspace))

# Example matrix
matrix = np.array([[1, 2, 3, 4],
                   [5, 6, 7, 8],
                   [9, 10, 11, 12],
                   [13, 14, 15, 16]])

# Shift all columns by one column using numpy.roll
shifted_matrix = np.roll(matrix, shift=-1, axis=1)

print("Original Matrix:")
print(matrix)
print("\nShifted Matrix:")
print(shifted_matrix)


# Plot the vector
plt.figure(figsize=(8, 6))
plt.plot(normalized_linspace, marker='o', linestyle='-')
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Plot of Vector')
plt.grid(True)
plt.show()

alphas = np.linspace(0,3,4)

print(matrix * alphas[:, np.newaxis].T)
print(np.sum(matrix * alphas[:, np.newaxis].T,axis=1))

print(np.append([0], np.logspace(-3,-0,10)))
print(np.logspace(0.0,1.0,11)/10)