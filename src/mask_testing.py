import numpy as np

# Données fictives
T = np.array([[1, 2],
              [3, 4],
              [5, 6]])

infected = np.array([1, 0])
memory_persuasion = np.array([0, 1])
infection_step = np.array([0, 1, 0])
i = 1

# Initialisation de la matrice mémoire
memory_matrix = np.zeros((3, 2))

# Application du masque
mask = infection_step == i - 1  # ici i-1 = 0, donc mask = [True, False, True]

# Calcul et mise à jour
memory_matrix[mask, -1] = T[mask].dot(infected + memory_persuasion)

print(memory_matrix)
