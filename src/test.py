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

# Plot the vector
plt.figure(figsize=(8, 6))
plt.plot(normalized_linspace, marker='o', linestyle='-')
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Plot of Vector')
plt.grid(True)
plt.show()