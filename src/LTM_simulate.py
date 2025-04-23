from utils import *
from network_functions import *

'''

'''

#* Paramètres de simulation
network_root = 'LTM'
network_class = ['ws']
N = [1000]
K = [8,16]
Tau = [0,1,2,4,8,16,32,64]
Persuasion = [0,0.01,0.05,0.1,0.2,0.5,1]
cascades = np.round(np.linspace(0.1,0.9,9),1)
threshold = np.linspace(0.01,0.5,16)
# Choix random de 5% pour le seed_nodes si percent_of_seeds = 0
# Si percent_of_seeds > 0, alors on prend le pourcentage de noeuds avec le plus grand degré
# Si percent_of_seeds < 0, alors on prend un pourcentage de noeuds avec le plus petit degré
percent_of_seeds = 0.05 # pourcentage de noeuds pour seed_nodes (exemple : int(n/20))

# Paramètres de dataframe
cols=['ID', 'network', 'p','th', 'seed']+ cascades.astype('str').tolist()
b = {'ID':[],'network':[], 'CC':[],'T':[],'p':[],'SP':[]}

for network_type in network_class:
    for n in N:        
        for k in K:
            network_props = pd.DataFrame(data=b)
            network_props.set_index(['ID','p'],inplace=True)
            df = pd.DataFrame(columns = cols)
                
            for t in Tau:
                for pers in Persuasion:
                    print()
                    print("Running simulation : network_type =",network_type,", N =",n,", K =",k,", t =",t,", pers =",pers)
                    print()
                    # Génère un vecteur de poids pour le mécanisme d'inertie
                    weights = gen_weights(t)
                    # The polirization file counter
                    count = 0
                    nets_prop_file =  f'{network_root}/{network_type}/{n}/{k}/{t}/{pers}/props.csv'
                    polarization_file = f'{network_root}/{network_type}/{n}/{k}/{t}/{pers}/polarization.csv'
                    ## If the network_path does not exist it is created
                    if not os.path.exists(f'{network_root}/{network_type}/{n}/{k}/{t}/{pers}'):
                        os.makedirs(f'{network_root}/{network_type}/{n}/{k}/{t}/{pers}')
                    
                    '''    
                    config_path = f'{network_root}/{network_type}/{n}/{k}/{t}/{pers}/config.ini'

                    config = configparser.ConfigParser()

                    config['Configuration'] = {'cascade': cascades, 'threshold' : threshold, 'probabilities' : probabilities, 'Weights' : weights, 'Nodes' : n, 'Neighbors' : k, 'Memory' : t, 'Network type' : network_type, 'Persuasion' : pers}

                    # Ecrire la configuration dans un fichier
                    with open(config_path, 'w') as configfile:
                        config.write(configfile)
                    '''
                        
                    #### Charger les graphs et lancer les simulations ####
                    for graph_path in list(get_recursive_graph_paths(f'Networks/{network_type}/{n}/{k}'))[:]:
                        start = time.perf_counter() # start time counter
                        G = gt.load_graph(str(graph_path))

                        G = get_local_clutsering(G)
                        G = get_transitivity(G)
                        G = get_ave_shortest_path(G)
                        G = get_laplacian_eigenvalues(G)
                        G = get_kirchhoff_index(G)
                        
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'network'] = 'ws'
                        # network_props.loc[(G.gp['ID'],G.gp['probability']),'network'] = 'mhk'
                        
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'CC'] = sum((G.vp.local_clustering.get_array()))/len(G.get_vertices())
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'T'] = G.gp.transitivity 
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'SP'] = G.gp.get('shortest_path')  
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'l2'] = np.sort(G.vp.eig_laplacian.a)[1]
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'lmax_l2'] =  np.max(G.vp.eig_laplacian.a) / np.sort(G.vp.eig_laplacian.a)[1]
                        network_props.loc[(G.gp['ID'],G.gp['probability']),'Rg'] =  n*np.sum(1/np.sort(G.vp.eig_laplacian.a)[1:])

                        #* Seed nodes selection
                        
                        #+ Harmonic centrality
                        # degrees = G.degree_property_map("total")
                        # max_degree_vertex = degrees.a.argmax()
                        # seed_nodes = np.array([max_degree_vertex])
                        ##harmonic_centrality = nx.harmonic_centrality(nx.from_numpy_array(gt.spectral.adjacency(G).T.toarray()))
                        # max_harmonic_node = max(harmonic_centrality, key=harmonic_centrality.get)
                        # seed_nodes = np.array([max_harmonic_node])
                        ##seed_nodes = sorted(harmonic_centrality, key=harmonic_centrality.get, reverse=True)[:10] # run que sur les 10 premiers noeuds avec la centralité la plus grande
                        # print(f'hc:{max_harmonic_node} VS nd:{max_degree_vertex}')
                        
                        #+ Degree centrality
                        degrees = dict(nx.from_numpy_array(gt.spectral.adjacency(G).T.toarray()).degree())
                        percent = int(percent_of_seeds * n)
                        if percent > 0:
                            sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)
                            seed_nodes = sorted_nodes[:percent]
                        elif percent < 0: #! erreur ? car selected_nodes[percent:] donne les 95% 
                            sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)
                            seed_nodes = sorted_nodes[percent:]
                        else:
                            #randomly select 5% of the nodes if you put '0' at percentage
                            all_nodes = list(degrees.keys())
                            num_nodes_to_choose = int(len(all_nodes) * 5 / 100)
                            seed_nodes = random.sample(all_nodes, num_nodes_to_choose)
                        #print("degrees : ", degrees.values())

                        for seed in seed_nodes:
                            #* Modèle complet (inertie + persuasion)
                            infected_vectormap, _, _ = linear_threshold_memory_model(G,threshold,pers,t,weights,seed_nodes=[seed])
                            spread = gt.ungroup_vector_property(infected_vectormap,range(len(threshold)))
                            data = np.empty((len(threshold),len(cascades),)) * np.nan
                            for idx,th in enumerate(threshold):
                                speeds = []
                                # speeds = np.empty((1,len(cascades),)) * np.nan
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
                                # print('th:',th)
                                nan_padding = len(cascades) - len(speeds)
                                speeds = np.pad(speeds, (0, nan_padding), constant_values=np.nan)
                                # print(speeds)
                                data[idx,:] = speeds
                                df.loc[count] = [G.gp.ID] + [G.gp.ntype] + [G.gp.probability] + [th] + [seed] + list(speeds)
                                count += 1
                        end = time.perf_counter()
                        print(f"Execution time for loading + running graph {G.gp.ID} : {end - start:.4f} secondes")
                                
                    df.to_csv(polarization_file, sep='\t',index = False)
                    network_props.to_csv(nets_prop_file,sep='\t',mode='w',header=True)