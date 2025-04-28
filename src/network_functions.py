from utils import *

def ws_network(N,k,p,seed=None):
    """
    Function for creating a Modified Holme-Kim network
    Takes inputs:
       N: int, Number of nodes
       k: integer, The mean degree of nodes
       p: Probability of rewireing of edges on the graph
       seed: Integer for the numpy random seed function
    Returnd:
       G: a graphtool graph

    """ 
    if seed is None:
        seed=np.random.randint(2**63)

    G = gt.Graph(directed=False)
    G.add_edge_list(np.transpose(sp.sparse.tril(nx.adjacency_matrix(nx.connected_watts_strogatz_graph(N,k,p,tries=1000000,seed=seed))).nonzero()))

    return G

def mhk_network(N,k,p,seed=None):
    """
    Function for creating a Modified Holme-Kim network
    Takes inputs:
       N: int, Number of nodes
       k: integer, The mean degree of nodes
       p: Probability of rewireing of edges on the graph
       seed: Integer for the numpy random seed function
    Returnd:
       G: a graphtool graph

    """ 
    if seed is None:
        seed=np.random.randint(2**63)
    G = nx.connected_watts_strogatz_graph(k,2,0)

    
    for n in range(k,N):
        anchor = np.random.choice(list(G.nodes()))   
        anchor_neigh = list(G.neighbors(anchor))

        for i in range(int(k)):

            if np.random.random() < p:     
                G.add_edge(np.random.choice(anchor_neigh),n)

            else: 
                try:
                    temp = np.random.choice(np.setdiff1d(G.nodes(),np.append(anchor_neigh,[anchor, n])))
                    G.add_edge(temp,n)
                except:
                    temp = np.random.choice(anchor_neigh)
                    G.add_edge(temp,n)
        G.add_edge(anchor,n)
    
    T = gt.Graph(directed=False)
    T.add_edge_list(np.transpose(sp.sparse.tril(nx.adjacency_matrix(G)).nonzero()))

    return T

def ke_network(n,m):

    #https://rf.mokslasplius.lt/acchieving-high-clustering-in-scale-free-networks/
    G = nx.connected_watts_strogatz_graph(m,m,0)
    active_nodes = list(G.nodes())
    for i in range(m,n):
        for k in active_nodes:
            G.add_edge(k,i)
        active_nodes.append(i)
        active_nodes.remove(random.choice(active_nodes))
        
    T = gt.Graph(directed=False)
    T.add_edge_list(np.transpose(sp.sparse.tril(nx.adjacency_matrix(G)).nonzero()))
    
#     w = np.logspace(-4,1,100)
#     network_type = 'ke'
#     network_path = 'ke_nets'
    
#     T.vertex_properties['gains'] = get_gain(T,w,n)
#     T.graph_properties['ID'] = T.new_graph_property('int64_t',val=int(m))
#     T.graph_properties['k'] = T.new_graph_property('int64_t',val=int(m))
#     T.graph_properties['ntype'] = T.new_graph_property('string',val=network_type)
#     frequencies = T.new_graph_property('vector<double>',val=w)
#     T.graph_properties['frequencies'] = frequencies
#     T = get_local_clutsering(T)
#     T = get_transitivity(T)
#     T = get_ave_shortest_path(T)
#     T = get_laplacian_eigenvalues(T)
    
    
    return T

@njit
def compute(memory_matrix,T,infected,memory_persuasion,weights,th,infection_step,i):
    #* Implémentation optimisée
    # Ne regarder que les infection_step qui ont la valeur i-1. Ceci signifie que les noeuds dont la valeur est i-1 à l'indice, viennent d'être infectés à l'itération précédente
    # Des noeuds fraichement infectés on veut regarder si ils parviennent à infecter d'autres noeuds. Il faut tout de même regarder les voisins de ces noeuds aussi...
    # Sans recalculer tous les noeuds !!!
    # Ne regarder que les voisins non-infectés des noeuds infectés
    #print("computing index_noeud_hors_infection")
    #index_noeud_hors_infection = ~np.any(T[infected == 1] != 0, axis=0) # indices des noeuds non-infectés pas en contact avec des noeuds infectés
    index_noeud_hors_infection = np.ones(len(infected), dtype=bool_)
    t_size = T[infected == 1].shape
    for colonne in range(t_size[1]):
        for ligne in range(t_size[0]):
            if T[infected == 1][ligne,colonne] != 0:
                index_noeud_hors_infection[colonne] = False 
                break  
    # np.full_like(infected, True, dtype=bool), init tout les noeuds à True pour le calcul
    #print("computing index_noeud_non_infectes_en_contact")
    index_noeud_non_infectes_en_contact = np.logical_xor(np.ones(len(infected), dtype=bool_),np.logical_or(infected.astype(bool_),index_noeud_hors_infection)) # noeuds non-infectés en contact avec des noeuds infectés
    #np.place(memory_matrix[:,-1],index_noeud_non_infectes_en_contact, T[index_noeud_non_infectes_en_contact].dot(infected + memory_persuasion))
    #print("computing memory_matrix")
    memory_matrix[index_noeud_non_infectes_en_contact,-1] = T[index_noeud_non_infectes_en_contact].dot(infected + memory_persuasion)
    #print("computing infected + infected_step")
    infected[memory_matrix @ weights >= th] = 1
    infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = i

def linear_threshold_memory_model(G,threshold,persuasion_step,tau,weights,seed_nodes=None,init_spread=True,max_iter=None):
    
    '''
    Mémoire de persuasion : Un vecteur de taille n est initialement rempli de 0.
    A chaque itération de la simulation les valeurs de tous les indexes dont les noeuds sont infectés, 
    sont incrémentés d'un pas défini par l'incrément "persuasion_step"
    
    Mémoire d'inertie : La variable tau représente le nombre d'itérations passées qui sont prises en compte.
    Initialement nous créons une matrice de n par tau qui équivaut à tau vecteurs avec chaque vecteur décrivant
    l'état d'infection de tous les noeuds du système. Initialement la matrice est remplie du vecteur de configuration 
    initial et à chaque itération les colonnes se décalent afin d'intégrer les nouveaux états calculés. Afin de calculer 
    l'état t+1 à partir de t nous prenons en compte tau vecteur décrivant les tau itérations dans le passé. Une pondération 
    déterminée par le vecteur "weights" permet la combinaison des différents vecteurs en un seul résultant.
    '''
    
    if seed_nodes == None:
        [seed_nodes for x in np.random.choice(G.get_vertices(),1)]

    if not type(seed_nodes) is list:
        seed_nodes = np.random.choice(G.get_vertices(),seed_nodes)

    if max_iter is None:
        max_iter = G.num_vertices()

    if not type(threshold) is list:
        [threshold]

    infections = []
    degree_dist = G.get_out_degrees(G.get_vertices())

    # T est la matrice d'adjacence normalisée par le degré des noeuds
    T = np.array((graph_tool.spectral.adjacency(G).T.toarray() / degree_dist).T)
    
    '''
    with open("matrice_T.txt", 'w') as f:
        for row in T:
            # Formater chaque valeur : 4 caractères, 1 chiffre après la virgule
            line = ' '.join(f"{val:4.2f}" for val in row)
            f.write(line + '\n')
            
    with open("matrice_T.T.txt", 'w') as f:
        for row in T.T:
            line = ' '.join(f"{val:4.2f}" for val in row)
            f.write(line + '\n')
    '''
    
    for th in threshold:
        infected = np.zeros(G.num_vertices(),dtype=int)
        infection_step = np.full(G.num_vertices(),np.inf,dtype=float)
        #node_list = np.arange(G.num_vertices(),dtype=int)

        #Infect the seed nodes
        infected[seed_nodes] = 1
        #Record seed nodes infected at t=-1
        infection_step[seed_nodes] = -1
        
        # Vecteur de taille n défini à 0
        memory_persuasion = np.zeros(G.num_vertices(),dtype=float)
        # Matrix of memory + current (last column)
        memory_matrix = np.zeros((G.num_vertices(),tau+1),dtype=float)

        '''
        T.dot(infected) : donne un vecteur de taille n avec la somme des états d'infection des voisins
        T.dot(infected) > 0 : renvoie un vecteur (booléen) de taille n avec 1 si le noeud est infecté et 0 sinon
        infected[T.dot(infected) > 0] = 1 : infecte les noeuds qui ont au moins un voisin infecté
        infected[T.dot(infected) >= th] = 1 : infecte les noeuds qui ont au moins th % de voisins infectés
        infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = i : définit à i (le pas de temps) les indices des noeuds infectés à l'itération i. Seul les valeurs à -inf ne peuvent être modifées
        '''

        #Initial spread, if choosen
        if init_spread:
            infected[T.dot(infected) > 0] = 1
            for i in range(tau+1):
                # Memory filled with initial configuration
                memory_matrix[:,i] = T.dot(infected)
            infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = 0
            i = 1
        else:
            i = 0
        '''
        Continue la simulation tant que : 
            1. tous les noeuds ne sont pas infectés
            2. le nombre d'itérations max n'est pas dépassé
            3. on n'est pas dans une situation de stagnation (pas d'infection à l'itération précédente)
        '''
        while (not np.all(infected) and (i < max_iter) and i-1 in infection_step):
            #print("iteration : ",i)
            #* Mécanisme de persuasion
            if(persuasion_step != 0): # la valeur de 0 n'active pas le mécanisme de mémoire persuasive
                memory_persuasion[infected != 0] += persuasion_step # incrémente si le noeuds à l'incide est infecté
            
            #* Mécanisme d'inertie
            if(tau != 0):
                # Shift columns one position left
                memory_matrix = np.roll(memory_matrix, shift=-1, axis=1)
                # Update last column of matrix with current state
                
            Opti = False
            if(not Opti):
                #* Standard implementation
                memory_matrix[:,-1] = T.dot(infected + memory_persuasion)
                test_val = memory_matrix @ weights
                infected[memory_matrix @ weights >= th] = 1
                infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = i
            else:
                #* Implémentation optimisée
                # Ne regarder que les infection_step qui ont la valeur i-1. Ceci signifie que les noeuds dont la valeur est i-1 à l'indice, viennent d'être infectés à l'itération précédente
                # Des noeuds fraichement infectés on veut regarder si ils parviennent à infecter d'autres noeuds. Il faut tout de même regarder les voisins de ces noeuds aussi...
                # Sans recalculer tous les noeuds !!!
                # Ne regarder que les voisins non-infectés des noeuds infectés
                #print("computing index_noeud_hors_infection")
                index_noeud_hors_infection = ~np.any(T[infected == 1] != 0, axis=0) # indices des noeuds non-infectés pas en contact avec des noeuds infectés
                # np.full_like(infected, True, dtype=bool), init tout les noeuds à True pour le calcul
                #print("computing index_noeud_non_infectes_en_contact")
                index_noeud_non_infectes_en_contact = np.logical_xor(np.full_like(infected, True, dtype=bool),np.logical_or(infected.astype(bool),index_noeud_hors_infection)) # noeuds non-infectés en contact avec des noeuds infectés
                #np.place(memory_matrix[:,-1],index_noeud_non_infectes_en_contact, T[index_noeud_non_infectes_en_contact].dot(infected + memory_persuasion))
                #print("computing memory_matrix")
                memory_matrix[index_noeud_non_infectes_en_contact,-1] = T[index_noeud_non_infectes_en_contact].dot(infected + memory_persuasion)
                #print("computing infected + infected_step")
                infected[memory_matrix @ weights >= th] = 1
                infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = i
                #compute(memory_matrix,T,infected,memory_persuasion,weights,th,infection_step,i)
            i += 1   
        infected_step = G.new_vp(value_type='int',vals=infection_step)
        infections.append(infected_step)

    infected_vectormap = gt.group_vector_property(infections)
    threshold_vector = G.new_gp(value_type='vector<double>',val=threshold)

    return infected_vectormap, seed_nodes, threshold_vector

def linear_threshold_model(G,threshold,seed_nodes=None,init_spread=True,max_iter=None):
    """
    Runs the linear threshold model on a grap_tool graph G. S_i(t+1) = {1 if <s_j>i~j > T otherwise 0. If the average state of neighbours of i is more than the threshold T, switch the state from 0 to 1

    Takes inputs:
       G: Graph_tool graph
       threshold: float or list of the thresholds as a fraction of neighbors that need to be active to transmit
       seed_nodes: int or list of nodes to start infection from, if None a single random agent is choosen. If int, a number of random nodes is choosen, if a list of nodes those are the initial seeds.
       init_spread: Bool, Wether to spread the infection to the neighbors of seed nodes
       max_iter = maximum number of iterations before stopping. If all nodes in the graph are active the model will stop
    Returnd:
       infected_step: Graph_tool VertexPropertyMap with an interger denoting the step of activation for a given vertex, None if never activated. Velues are stored as vectors, with each index corresponding to a given threshold.
       seed_nodes: Graph_tool GraphProperty with the initail seed nodes for the model run

    """

    if seed_nodes == None:
        [seed_nodes for x in np.random.choice(G.get_vertices(),1)]

    if not type(seed_nodes) is list:
        seed_nodes = np.random.choice(G.get_vertices(),seed_nodes)

    if max_iter is None:
        max_iter = G.num_vertices()

    if not type(threshold) is list:
        [threshold]

    infections = []
    degree_dist = G.get_out_degrees(G.get_vertices())

    T = np.array((graph_tool.spectral.adjacency(G).T.toarray() / degree_dist).T)  

    for th in threshold:
        # Choose the initial infected nodes
        infected = np.zeros(G.num_vertices(),dtype=int)
        infection_step = np.full(G.num_vertices(),np.inf,dtype=float)
        node_list = np.arange(G.num_vertices(),dtype=int)

        #Infect the seed nodes
        infected[seed_nodes] = 1
        #Record seed nodes infected at t=-1
        infection_step[seed_nodes] = -1

        #Initial spread, if choosen
        if init_spread:
            infected[T.dot(infected) > 0] = 1
            infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = 0
            i = 1
        else:
            i = 0
        while (not all(infected) and (i < max_iter) and i-1 in infection_step):
            infected[T.dot(infected) >= th] = 1
            infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = i
            i += 1
        infected_step = G.new_vp(value_type='int',vals=infection_step)
        infections.append(infected_step)

    infected_vectormap = gt.group_vector_property(infections)
    threshold_vector = G.new_gp(value_type='vector<double>',val=threshold)
    # G.vp['infected_step'] = infected_vectormap    
    # G.gp['threshold_vector'] = G.new_gp(value_type='vector<double>',val=threshold)
    # G.gp['cascades'] = G.new_gp(value_type='vector<double>',val=cascades)
    # G.gp['seed_nodes'] = G.new_gp(value_type='vector<double>',val=[seed_nodes])
    
    return infected_vectormap, seed_nodes, threshold_vector  

def gen_weights(tau):
    # Weight vector
    # Define weights values (different methods)
    # 1 : linear
    weights = np.linspace(1/(tau+1), 1, tau+1)
    # 2 : log
    #weights = np.logspace(1/(tau+1), 1, tau+1, base=10)
    # 3 : uniform
    #weights = np.full(tau+1,1)
    # Reverse alpha vector
    #weights = np.flip(weights)
    
    #? Test : weight for current state = 1 the rest = 0
    #weights = np.full(tau+1,0)
    #weights[-1] = 1 # last alpha = current state
    
    # Normalization
    weights = weights / np.sum(weights)
    #print("sum of weights is : ",np.sum(weights))
    return weights

def get_gain(G,graph,w,N):
    L = gt.spectral.laplacian(graph,norm=False) #the build in normalized gives the symetric normalized laplacian, but we want the random walk normalized laplacian
    
    L = (L/L.diagonal()).T  ## Random walk normalization  D^-1 L = LD^-1 because L is symetric

    L = L.toarray()
    h2 = G.new_vertex_property('vector<double>')
    for g in G.vertices():
        ida = np.arange(N) != g
        idb = np.arange(N) == g
        A = L[np.ix_(ida,ida)].astype(complex)
        B = L[np.ix_(ida,idb)]

        H2 = []    
        for f in w:
            np.fill_diagonal(A,1.0+1j*f)
            #A.setdiag(f*1j-1)
            h = linalg.solve(A,-B)
            
            H2.append(linalg.norm(h)**2)

        h2[g] = H2

    return h2

def get_laplacian_eigenvalues(G):
    """
    """

    if not G.vertex_properties.get('eig_laplacian',False):
        eig_lap = np.linalg.eigvalsh(gt.spectral.laplacian(G,norm=False).todense())
        G.vp['eig_laplacian'] =  G.new_vertex_property('double',vals=eig_lap)
    return G

def get_kirchhoff_index(G):
    G = get_laplacian_eigenvalues(G)
    G.graph_properties['kirchhoff'] = G.new_graph_property('int64_t',sum(1/np.sort(G.vp.eig_laplacian.get_array())[1:]))
    return G

def get_recursive_graph_paths(root):
    paths = Path(Path(root)).rglob('*.gt')
    return paths

def get_local_clutsering(G):
    if not G.vertex_properties.get('local_clustering',False):
        G.vertex_properties['local_clustering'] = graph_tool.clustering.local_clustering(G)
    return G

def get_transitivity(G):
    if not G.gp.get('transitivity',False):
        trans = G.new_graph_property('double',val=graph_tool.clustering.global_clustering(G)[0])
        G.graph_properties['transitivity'] =  trans
    return G

def get_ave_shortest_path(G):
    if not G.gp.get('shortest_path',False):
        G.gp['shortest_path'] =  G.new_graph_property('double',val=np.sum( graph_tool.topology.shortest_distance(G).get_2d_array(range(G.num_vertices())))/(G.num_vertices()*(G.num_vertices()-1)))
    return G