import numpy as np
import networkx as nx
import scipy as sp
from pathlib import Path
from scipy import sparse
import graph_tool as gt
import graph_tool.clustering
import graph_tool.spectral
import graph_tool.topology
from pathlib import Path
from scipy import linalg
import random


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

def linear_threshold_memory_model(G,threshold,tau,alphas,seed_nodes=None,init_spread=True,max_iter=None):

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
    #print(graph_tool.spectral.adjacency(G).T.toarray().shape) 
    #print(degree_dist.shape)

    for th in threshold:
        # Matrix of memory + current (last column)
        memory_matrix = np.zeros((G.num_vertices(),tau+1),dtype=float)
        # Choose the initial infected nodes
        infected = np.zeros(G.num_vertices(),dtype=int)
        #memory = np.zeros(G.num_vertices(),dtype=float)
        #current = np.zeros(G.num_vertices(),dtype=float)
        infection_step = np.full(G.num_vertices(),np.inf,dtype=float)
        node_list = np.arange(G.num_vertices(),dtype=int)

        #Infect the seed nodes
        infected[seed_nodes] = 1
        #Record seed nodes infected at t=-1
        infection_step[seed_nodes] = -1

        #Initial spread, if choosen
        if init_spread:
            # Treshold defined to 0 for init step
            infected[T.dot(infected) > 0] = 1
            for i in range(tau+1):
                # Memory filled with initial configuration
                memory_matrix[:,i] = T.dot(infected)
            infection_step[np.logical_and(infected > 0, np.isinf(infection_step))] = 0
            i = 1
        else:
            i = 0
        while (not all(infected) and (i < max_iter) and i-1 in infection_step):
            # Memorise <Sj> of iteration t-1
            infected[(np.sum(memory_matrix * alphas[:, np.newaxis].T,axis=1)) >= th] = 1
            # Shift columns one position left
            memory_matrix = np.roll(memory_matrix, shift=-1, axis=1)
            # Update last column of matrix with current state
            memory_matrix[:,-1] = T.dot(infected)
            #? Test : Update all columns with current state
            #for j in range(tau+1):
            #    memory_matrix[:,j] = T.dot(infected)
            
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

def gen_alphas(tau):
    # Alphas vector
    # Define alphas values (different methods)
    # 1 : linear
    alphas = np.linspace(1/(tau+1), 1, tau+1)
    # 2 : log
    #alphas = np.logspace(1/(tau+1), 1, tau+1, base=10)
    # 3 : uniform
    #alphas = np.full(tau+1,1)
    # Reverse alpha vector
    #alphas = np.flip(alphas)
    
    #? Test : alpha for current state = 1 the rest = 0
    #alphas = np.full(tau+1,0)
    #alphas[-1] = 1 # last alpha = current state
    
    # Normalization
    alphas = alphas / np.sum(alphas)
    #print("sum of alphas is : ",np.sum(alphas))
    return alphas


def get_gain(graph,w,N):
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