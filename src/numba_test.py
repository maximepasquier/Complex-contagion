from numba import njit  # njit = no python mode = maximum perf

@njit
def somme_cond(liste, theta):
    n = len(liste)
    res = [0] * n
    for i in range(n):
        if liste[i] > theta:
            res[i] = 1
    return res

# Exemple d'utilisation
import numpy as np
liste = np.random.rand(1000000)
theta = 0.5

out = somme_cond(liste, theta)
print(out[:10])  # Affiche les 10 premiers éléments du résultat