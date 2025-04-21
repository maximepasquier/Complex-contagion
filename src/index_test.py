import numpy as np

# Exemple de matrice
A = np.array([
    [0, 1, 0],
    [0, 0, 0],
    [0, 2, 0]
])

# Teste si toutes les valeurs d'une colonne sont 0
zero_columns_mask = np.all(A == 0, axis=0)

print(zero_columns_mask)
