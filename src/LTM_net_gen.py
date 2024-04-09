'''
Generates and Saves graphs in gt format and generates 2 csv files (props and polarization)
'''

import numpy as np
import networkx as nx
import scipy as sp
from pathlib import Path
from scipy import sparse
import graph_tool.clustering
import graph_tool.spectral
import graph_tool.topology
import matplotlib.pyplot as plt
import graph_tool as gt
import glob
import pandas as pd
import os
from scipy import linalg
from scipy.sparse import coo_array
import time
from matplotlib.colors import to_hex
import matplotlib
import sys
import warnings
import random
warnings.filterwarnings("ignore")
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble = r'\usepackage{mathptmx}')

from network_functions import *

network_root = 'LTM_networks'
network_type = 'mhk'
N = [1000]
K = [16]

cascades = np.round(np.linspace(0.1,0.9,9),1)
threshold = np.linspace(0.01,0.5,16)    


probabilities = np.append([0], np.linspace(0.9,1,10))
# probabilities = [0]
cols=['ID', 'network', 'p','th', 'seed']+ cascades.astype('str').tolist()
b = {'ID':[],'network':[], 'CC':[],'T':[],'p':[],'SP':[]}

for n in N:
    for k in K:
        
        # The polirization file counter
        count = 0
        nets_prop_file =  f'{network_root}/{network_type}/{n}/{k}_props.csv'
        polarization_file = f'{network_root}/{network_type}/{n}/polarization_{k}.csv'
        network_path = f"{network_root}/{network_type}/{n}/{k}"
                ## If the network_path does not exist it is created
        if not os.path.exists(network_path):
            os.makedirs(network_path)

        network_props = pd.DataFrame(data=b)
        network_props.set_index(['ID','p'],inplace=True)
        df = pd.DataFrame(columns = cols)

        for prob in probabilities:

            G = mhk_network(n,k,prob)


            G.graph_properties['ID'] = G.new_graph_property('int64_t',val=int(time.time()*1000))
            G.graph_properties['ntype'] = G.new_graph_property('string',val=network_type)
            G.graph_properties['probability'] = G.new_graph_property('double',prob)
            G.graph_properties['cascades'] = G.new_gp(value_type='vector<double>',val=cascades)
            G = get_local_clutsering(G)
            G = get_transitivity(G)
            G = get_ave_shortest_path(G)
            G = get_laplacian_eigenvalues(G)
            G = get_kirchhoff_index(G)
            G.save(f'{network_path}/{G.gp.ID}.gt')
            
            network_props.loc[(G.gp['ID'],G.gp['probability']),'network'] = 'mhk'
            network_props.loc[(G.gp['ID'],G.gp['probability']),'CC'] = sum((G.vp.local_clustering.get_array()))/len(G.get_vertices())
            network_props.loc[(G.gp['ID'],G.gp['probability']),'T'] = G.gp.transitivity 
            network_props.loc[(G.gp['ID'],G.gp['probability']),'SP'] = G.gp.get('shortest_path')  
            network_props.loc[(G.gp['ID'],G.gp['probability']),'l2'] = np.sort(G.vp.eig_laplacian.a)[1]
            network_props.loc[(G.gp['ID'],G.gp['probability']),'lmax_l2'] =  np.max(G.vp.eig_laplacian.a) / np.sort(G.vp.eig_laplacian.a)[1]
            network_props.loc[(G.gp['ID'],G.gp['probability']),'Rg'] =  n*np.sum(1/np.sort(G.vp.eig_laplacian.a)[1:])


            # degrees = G.degree_property_map("total")
            # max_degree_vertex = degrees.a.argmax()
            # seed_nodes = np.array([max_degree_vertex])
            harmonic_centrality = nx.harmonic_centrality(nx.from_numpy_array(gt.spectral.adjacency(G).T.toarray()))
            # max_harmonic_node = max(harmonic_centrality, key=harmonic_centrality.get)
            # seed_nodes = np.array([max_harmonic_node])
            seed_nodes = sorted(harmonic_centrality, key=harmonic_centrality.get, reverse=True)[:100]
            # print(f'hc:{max_harmonic_node} VS nd:{max_degree_vertex}')

            print('G.ID:',G.gp.ID)

            for seed in seed_nodes:
                infected_vectormap, selected_seed, _ = linear_threshold_model(G,threshold,seed_nodes=[seed])

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
        df.to_csv(polarization_file, sep='\t',index = False)
        network_props.to_csv(nets_prop_file,sep='\t',mode='w',header=True)