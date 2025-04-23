from utils import *
from network_functions import *

'''
Génère les réseaux en fonction des paramètres
    - network_root : nom du dossier contenant les réseaux
    - network_class : liste des classes de réseau
    - N : liste de nombre de noeuds
    - K : liste de nombre de voisins moyen
    - probabilities : liste de probabilités de rewiring
Crée uniquement les réseaux absent dans les sous-dossiers
'''

#* Paramètres de réseau
network_root = 'LTM'
probabilities = np.append([0], np.logspace(-3,-0,10))
probabilities = probabilities[[0,5,7,8]] # filtre sur les valeurs de probabilités
network_class = ['ws']
N = [1000]
K = [8,16]

#* Création des réseaux
for network_type in network_class:
    for n in N:        
        for k in K:
            # Crée le dossier pour les réseaux
            if not os.path.exists(f'{network_root}/{network_type}/Networks/{n}/{k}'):
                os.makedirs(f'{network_root}/{network_type}/Networks/{n}/{k}')
            # Make a dict with the probabilities as keys, to count the available networks
            graph_realized=dict.fromkeys(probabilities,0)
            # Count the graphs in network path and create missing networks according to  parameters (N,k,p)
            for graph_path in get_recursive_graph_paths(f'{network_root}/{network_type}/Networks/{n}/{k}'):
                g_old = gt.load_graph(str(graph_path))
                graph_realized[g_old.gp.probability] = 1
            #### Create missing networks ####
            for p in probabilities:
                if graph_realized[p] == 1: # déjà un graph créé pour ce rewiring
                    continue
                G = ws_network(n,k,p,seed=None)
                # Add relevant creation properties with the graph.
                # Nomme les graph en fonction de p
                G.graph_properties['ID'] = G.new_graph_property('string',val="Network_p="+str(p))
                G.graph_properties['ntype'] = G.new_graph_property('string',val=network_type)
                G.graph_properties['probability'] = G.new_graph_property('double',p)
                #G.graph_properties['cascades'] = G.new_gp(value_type='vector<double>',val=cascades)
                
                print("New graph created with ID: ", G.gp.ID)

                G.save(f'{network_root}/{network_type}/Networks/{n}/{k}/{G.gp.ID}.gt')