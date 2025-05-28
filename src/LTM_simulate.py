from utils import *
from network_functions import *

'''
Input :
    - network_root : nom du dossier contenant les réseaux
    - network_class : liste des classes de réseau
    - N : liste de nombre de noeuds
    - K : liste de nombre de voisins moyen
    - Persuasion : liste de biais de persuasion
    - Tau : liste de taille de la mémoire
    - cascades : liste de taille de cascade
    - threshold : liste de seuils
    - percent_of_seeds : pourcentage de noeuds pour seed_nodes
    
Output :
    - Arborescence de dossiers contenant les résultats des simulations (polarization.csv et props.csv)
    - Polarization.csv : fichier contenant les vitesses de polarisations
    - Props.csv : fichier contenant les propriétés du réseau

Le script run les analyses sur les réseaux défini par les paramètres de simulations ainsi que ceux générés par le script LTM_net_gen.py.
Des listes de N, K, Tau et Persuasion, toutes les cominaison de configurations sont générées et exécutées. 
Tous les résultats (polarization.csv et props.csv) sont stockés dans : LTM/{network_type}/{N}/{K}/{Tau}/{Persuasion}
'''  

#* Paramètres de simulation
network_root = 'LTM' # dossier de sortie
network_class = ['ws'] # type de réseaux
N = [500] # nombre de noeuds
K = [16] # nombre de voisins moyen
#waiting_counter_max = [0] # nombre max d'itérations admise entre deux pas d'évolution : 0 signifie consécutif
Tau = [0] # taille de la mémoire
cascades = np.round(np.linspace(0.1,0.9,9),1) # tailles des cascades
threshold = np.linspace(0.01,0.5,16) # thresholds
'''
Si top, alors on prend le pourcentage de noeuds avec le plus grand degré
Si bot, alors on prend le pourcentage de noeuds avec le plus petit degré
Si rand, alors on prend le pourcentage des noeuds au hasard
'''
seed_types = ['top','bot','rand'] # types de seed_nodes
fraction_of_seeds = 0.05 # fraction de noeuds pour seed_nodes
optimisation = False # True si on veut "optimiser" le code
memory_saturation = False # Si True alors la mémoire tau est saturée par l'état initial au début de la simulation

#* Paramètres de dataframe
cols=['ID', 'network', 'p','th', 'seed']+ cascades.astype('str').tolist()
b = {'ID':[],'network':[], 'CC':[],'T':[],'p':[],'SP':[]}

for network_type in network_class:
    for n in N:
        nb_seed_nodes = int(n * fraction_of_seeds) # nombre de seed_nodes        
        for k in K: 
            network_props = pd.DataFrame(data=b)
            network_props.set_index(['ID','p'],inplace=True)    
            for t in Tau:
                ## If the network_path does not exist it is created
                if not os.path.exists(f'{network_root}/{network_type}/{n}/{k}/{t}'):
                    os.makedirs(f'{network_root}/{network_type}/{n}/{k}/{t}')
                nets_prop_file =  f'{network_root}/{network_type}/{n}/{k}/{t}/props.csv'
                for seed_type in seed_types:
                    print()
                    print("Running simulation : network_type =",network_type,", N =",n,", K =",k,", t =",t,", seed_type =",seed_type)
                    df = pd.DataFrame(columns = cols) # dataframe pour les vitesses de polarisation
                    polarization_file = f'{network_root}/{network_type}/{n}/{k}/{t}/{seed_type}{fraction_of_seeds}polarization.csv'
                    # The polirization file counter
                    count = 0
                    #### Charger les graphs et lancer les simulations ####
                    for graph_path in list(get_recursive_graph_paths(f'Networks/{network_type}/{n}/{k}'))[:]:
                        start = time.perf_counter() # start time counter
                        G = gt.load_graph(str(graph_path))

                        G = get_local_clutsering(G)
                        G = get_transitivity(G)
                        G = get_ave_shortest_path(G)
                        G = get_laplacian_eigenvalues(G)
                        G = get_kirchhoff_index(G)
                        
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'network'] = network_type
                        
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'CC'] = sum((G.vp.local_clustering.get_array()))/len(G.get_vertices())
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'T'] = G.gp.transitivity 
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'SP'] = G.gp.get('shortest_path')  
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'l2'] = np.sort(G.vp.eig_laplacian.a)[1]
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'lmax_l2'] =  np.max(G.vp.eig_laplacian.a) / np.sort(G.vp.eig_laplacian.a)[1]
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'Rg'] =  n*np.sum(1/np.sort(G.vp.eig_laplacian.a)[1:])
                        
                        #+ Harmonic centrality
                        '''
                        degrees = G.degree_property_map("total")
                        max_degree_vertex = degrees.a.argmax()
                        seed_nodes = np.array([max_degree_vertex])
                        #harmonic_centrality = nx.harmonic_centrality(nx.from_numpy_array(gt.spectral.adjacency(G).T.toarray()))
                        max_harmonic_node = max(harmonic_centrality, key=harmonic_centrality.get)
                        seed_nodes = np.array([max_harmonic_node])
                        #seed_nodes = sorted(harmonic_centrality, key=harmonic_centrality.get, reverse=True)[:10] # run que sur les 10 premiers noeuds avec la centralité la plus grande
                        print(f'hc:{max_harmonic_node} VS nd:{max_degree_vertex}')
                        '''
                        #+ Degree centrality
                        degrees = dict(nx.from_numpy_array(gt.spectral.adjacency(G).T.toarray()).degree())
                        #print(degrees.values())
                        #* Seed nodes selection
                        if seed_type == 'top':
                            sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)
                            seed_nodes = sorted_nodes[:nb_seed_nodes]
                        elif seed_type == 'bot':
                            sorted_nodes = sorted(degrees, key=degrees.get, reverse=False)
                            seed_nodes = sorted_nodes[:nb_seed_nodes]
                        elif seed_type == 'rand':
                            all_nodes = list(degrees.keys())
                            seed_nodes = random.sample(all_nodes, nb_seed_nodes)

                        for seed in seed_nodes:
                            #* Modèle complet (inertie + persuasion)
                            infected_vectormap, _, _ = linear_threshold_memory_model(G,threshold,t,optimisation,memory_saturation,seed_nodes=[seed])
                            spread = gt.ungroup_vector_property(infected_vectormap,range(len(threshold)))
                            data = np.empty((len(threshold),len(cascades),)) * np.nan
                            for idx,th in enumerate(threshold):
                                speeds = [] # liste pour stocker les vitesses de polarisation
                                cascade_sizes = cascades ## Copy where we can shorten list on iterations in while loop
                                infected = 0             ## Reset on loop start

                                ## Find nodes infected at given step
                                val,counts = np.unique(spread[idx].a,return_counts=True)
                                ## Refactor to fraction of nodes
                                counts = counts / G.num_vertices()
                                ## For each step in LTM add newly infected to total
                                for i,new in enumerate(counts[val>-2]):
                                    infected += new
                                    ## Current number of infected are used to find polarization speed if larger than cascade size, and the step is after initialization, and all cascade sizes has not been exceeded yet.
                                    while len(cascade_sizes) > 0 and infected > cascade_sizes[0] and val[i] > 0:
                                        speeds.append(infected/val[i])
                                        cascade_sizes = cascade_sizes[1:]
                                nan_padding = len(cascades) - len(speeds)
                                speeds = np.pad(speeds, (0, nan_padding), constant_values=np.nan)
                                data[idx,:] = speeds
                                df.loc[count] = [G.gp.ID] + [G.gp.ntype] + [G.gp.probability] + [th] + [seed] + list(speeds)
                                count += 1
                        end = time.perf_counter()
                        print(f"Execution time for loading and running graph {G.gp.ID} : {end - start:.4f} secondes")
                        
                        '''
                        #* Select que la dernière seed (1 seule simulation) pour l'écriture dans les json
                        data_json = {}
                        spread = gt.ungroup_vector_property(infected_vectormap,range(len(threshold)))
                        for idx,th in enumerate(threshold):
                            val,counts = np.unique(spread[idx].a,return_counts=True)
                            data_json[str(th)] = {"val": val.tolist(), "counts": counts.tolist()}
                        waiting_dict_serializable = {
                            key: [arr.tolist() for arr in value] for key, value in waiting_dict.items()
                        }
                        #* Écrire dans un fichier JSON
                        with open(f'{network_root}/{network_type}/{n}/{k}/{t}/infected_p={G.gp.ID}.json', "w") as f:
                            json.dump(data_json, f, indent=2)
                        with open(f'{network_root}/{network_type}/{n}/{k}/{t}/waiting_dct_p={G.gp.ID}.json', 'w') as f:
                            json.dump(waiting_dict_serializable, f, indent=2)
                        '''
                            
                    df.to_csv(polarization_file, sep='\t',index = False)
                    network_props.to_csv(nets_prop_file,sep='\t',mode='w',header=True)