import numpy as np

A = np.array([
    [1,2],
    [4,6]
])  # shape (2, 3)

weights = np.array([1,2])  # shape (3,)

print(A @ weights)