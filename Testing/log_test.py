import numpy as np
import matplotlib.pyplot as plt

'''
Ce script génère et normalise des poids linéaires et logarithmiques pour une période donnée (tau).
Ca a pour but de visualiser les poids en fonction de tau, en utilisant des poids linéaires et logarithmiques.
'''

# Paramètres
tau = 8

# Poids linéaires
weights_lin = np.linspace(1 / (tau + 1), 1, tau + 1)
weights_lin = weights_lin / np.sum(weights_lin)

# Poids logarithmiques
weights_log = np.logspace(1 / (tau + 1), 1, tau + 1, base=20)
weights_log = weights_log / np.sum(weights_log)

print(weights_lin.sum())
print(weights_log.sum())

# Plot
plt.figure(figsize=(8, 4))
plt.plot(weights_lin, marker='o', label='Linéaire')
plt.plot(weights_log, marker='s', label='Logarithmique')
plt.title(f'Poids normalisés pour tau = {tau}')
plt.xlabel('Indice')
plt.ylabel('Poids')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
