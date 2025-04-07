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
network_root = 'LTM_networks'
#+ Watts_strogatz
probabilities = np.append([0], np.logspace(-3,-0,10)) # probabilité de rewireing
network_class = ['ws']
#+ Modified Holme-Kim
#network_class = ['mhk']

#* Paramètres de réseau
N = [100,200]
K = [4,8]
Tau = [0,2]
Persuasion = [0,0.05]

#* Paramètres de simulation
cascades = np.round(np.linspace(0.1,0.9,9),1)
threshold = np.linspace(0.01,0.5,16)

## How many realizations to do of each set of parameters
#desired_realizations= 1
## How many unique starting points to run the LTM from on a network
#unique_network_seeds = 1

#cols=['ID', 'network', 'p','th', 'seed']+ cascades.astype('str').tolist()
#b = {'ID':[],'network':[], 'CC':[],'T':[],'p':[],'SP':[]}

#* Simulations

ltm = LTM_memory(network_root,network_class,N,K,Tau,Persuasion,cascades,threshold,probabilities)

#ltm.visualize()
ltm.run()