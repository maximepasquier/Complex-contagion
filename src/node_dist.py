from collections import Counter
from utils import *
from network_functions import *

'''
Produit les bar plots de la distribution des degrés des noeuds des graphes.

Exemple avec :

mhk
N = 1000
K = 16
p = 0.5,0.95,0.99,1.0

Output :
figs/node_dist/...
'''

network_root = 'Networks'
network_class = ['mhk']
N = [1000]
K = [16]

#* Plotting
for network_type in network_class:               
    for n_nodes in N:
        for neighbor_k in K:
            #### Charger les graphs et lancer les simulations ####
            for graph_path in list(get_recursive_graph_paths(f'{network_root}/{network_type}/{n_nodes}/{neighbor_k}'))[:]:
                G = gt.load_graph(str(graph_path))
                G = get_local_clutsering(G)
                G = get_transitivity(G)
                G = get_ave_shortest_path(G)
                G = get_laplacian_eigenvalues(G)
                G = get_kirchhoff_index(G)
                
                #+ Degree centrality
                degrees = dict(nx.from_numpy_array(gt.spectral.adjacency(G).T.toarray()).degree()) 
                degree_values = list(degrees.values())

                # Binning des degrés
                bin_size = 5  # ajustable selon les cas
                max_degree = max(degree_values)
                binned_counts = Counter()
                
                for deg in degree_values:
                    bin_min = (deg // bin_size) * bin_size
                    bin_max = bin_min + bin_size - 1
                    bin_label = f"{bin_min}-{bin_max}"
                    binned_counts[bin_label] += 1

                # Trier les bins
                bin_labels = sorted(binned_counts.keys(), key=lambda s: int(s.split('-')[0]))
                bin_values = [binned_counts[label] for label in bin_labels]

                # Création du bar plot
                plt.bar(bin_labels, bin_values, color='skyblue', edgecolor='black')
                plt.yscale('log')
                plt.xlabel('Degré (regroupé en intervalles)')
                plt.ylabel('Nombre de noeuds')
                plt.title(f'Distribution des degrés \n Réseau {network_type} avec N={n_nodes}, K={neighbor_k}, clustering={G.vp.local_clustering.get_array().mean():.2f}')
                plt.xticks(bin_labels[::5])
                plt.grid(axis='y', linestyle='--', alpha=0.7)

                # Création du dossier de figures si besoin
                fig_path = f'figs/node_dist/{network_type}/{n_nodes}/{neighbor_k}'
                if not os.path.exists(fig_path):
                    os.makedirs(fig_path)
                
                # Sauvegarde de la figure
                plt.savefig(fig_path + f"/{G.gp['ID']}.pdf")
                plt.clf()
                print('fig saved')

