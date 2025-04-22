'''
Fichier main permettant de lancer les simulations LTM en déterminant les paramètres.

Input : 
    - Type de réseau
    - Nombre de noeuds
    - Nombre de voisins (moyens ?)
    - Nombre d'itérations pour processus markovien ou non
    - Terme de persuasion
    - Valeurs des cascades
    - Valeurs des seuils
    - Probabilité de "rewirering" pour les réseaux ws

Output : 
    - Génère les graphes en format gt et les sauvergarde
    - Génère les fichiers polarization.csv, props.csv et config.ini
    - Génère les figures pour l'ensemble des paramètres
'''

from utils import *
from network_functions import *
from LTM_memory_class import *


#* Type de réseau
network_root = 'LTM_networks' # nom du dossier contenant les réseaux + résultats des simulations
#+ Watts_strogatz
probabilities = np.append([0], np.logspace(-3,-0,10)) # vecteur de probabilités de rewiring (pour de multiples réseaux ws)
probabilities = probabilities[[0,5,7,8]] # filtre sur les valeurs de probabilités
#probabilities = probabilities[[5]]
network_class = ['ws'] # liste des classes de réseau
#+ Modified Holme-Kim
#network_class = ['mhk']

#* Paramètres de réseau
# Le nombre de configurations équivaut au produit de la longueur de l'ensemble des paramètres de réseau
N = [1000]          # nombre de noeuds
K = [8,16]            # nombre de voisins moyen
Tau = [0,1,2,4,8,16]           # nombre d'itérations prises en compte dans le mécanisme d'inertie
Persuasion = [0,0.001,0.01,0.05,0.1,0.2]    # valeurs d'influence pour la mécanisme de persuasion

#* Paramètres de simulation
cascades = np.round(np.linspace(0.1,0.9,9),1)   # liste de valeurs de cascades
threshold = np.linspace(0.01,0.5,16)            # liste de valeurs de seuil
#threshold = threshold[[5]]
# Choix random de 5% pour le seed_nodes si percent_of_seeds = 0
# Si percent_of_seeds > 0, alors on prend le pourcentage de noeuds avec le plus grand degré
# Si percent_of_seeds < 0, alors on prend un pourcentage de noeuds avec le plus petit degré
percent_of_seeds = 0.05 # pourcentage de noeuds pour seed_nodes (exemple : int(n/20))

#* Simulations
ltm = LTM_memory(network_root,network_class,N,K,Tau,Persuasion,cascades,threshold,probabilities,percent_of_seeds)
ltm.generate()
ltm.visualize()
ltm.analyse()