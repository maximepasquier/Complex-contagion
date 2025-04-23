import numpy as np
import matplotlib.pyplot as plt

def g(a,b):
    return a/b

def f(a,b,c):
    return (a+c)/(b+c)

a = np.linspace(0,99,100)
b = 100
c1 = np.linspace(0,9.9,100) # biais de persuasion = 0.1
c2 = np.linspace(0,19.9,100) # biais de persuasion = 0.2
c3 = np.linspace(0,29.9,100) # biais de persuasion = 0.3
c4 = np.linspace(0,39.9,100) # biais de persuasion = 0.4
c5 = np.linspace(0,49.9,100) # biais de persuasion = 0.5
c6 = np.linspace(0,59.9,100) # biais de persuasion = 0.6
c7 = np.linspace(0,69.9,100) # biais de persuasion = 0.7
c8 = np.linspace(0,79.9,100) # biais de persuasion = 0.8
c9 = np.linspace(0,89.9,100) # biais de persuasion = 0.9


# Création du plot
plt.figure(figsize=(10, 6))
plt.plot(a, g(a,b), label="base", color="red")
plt.plot(a, f(a, b, c1), label="biais = 0.1", color="tab:blue")
plt.plot(a, f(a, b, c2), label="biais = 0.2", color="tab:orange")
plt.plot(a, f(a, b, c3), label="biais = 0.3", color="tab:green")
plt.plot(a, f(a, b, c4), label="biais = 0.4", color="tab:cyan")
plt.plot(a, f(a, b, c5), label="biais = 0.5", color="tab:purple")
plt.plot(a, f(a, b, c6), label="biais = 0.6", color="tab:brown")
plt.plot(a, f(a, b, c7), label="biais = 0.7", color="tab:pink")
plt.plot(a, f(a, b, c8), label="biais = 0.8", color="tab:gray")
plt.plot(a, f(a, b, c9), label="biais = 0.9", color="tab:olive")
plt.title("Visualisation de l'influence des biases de persuasion sur le nombre d'infectés")
plt.xlabel("itérations")
plt.ylabel("% d'infectés")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
