#!/usr/bin/env python3
## Script for running the generating Small-world graphs, analyzing their network properties and running the Linear Threshold Model
import os
import numpy as np
import time
import graph_tool as gt
import random
from pathlib import Path

from network_functions import ws_network
from network_functions import get_recursive_graph_paths
from network_functions import linear_threshold_model
from network_functions import linear_threshold_memory_model
from network_functions import get_local_clutsering
from network_functions import get_transitivity
from network_functions import get_laplacian_eigenvalues
from network_functions import get_kirchhoff_index
from network_functions import get_ave_shortest_path
## File used to generate netwokrs, and subsequently do the full analysis
old_graph_list = []
#Parameters for networks
network_type = 'ws'
nodes=[1000]
average_degree=[16]
## Thresholds are chosen such that they approximately give the thresholds of [1,...,8]
thresholds = [0.3]

probabilities = [0.1]
## Root directory for storing networks
network_root = 'networks_test/'

##Cascade sizes of interest
cascades = np.round(np.linspace(0.1,0.9,9),1)

#Seed path for saving networks, put the
network_path = f"{network_root}{nodes[0]}/{average_degree[0]}"
## If the network_path does not exist it is created
if not os.path.exists(network_path):
    os.makedirs(network_path)

G = ws_network(nodes[0],average_degree[0],probabilities[0],seed=None)
# Add relevant creation properties with the graph.
# Name graphs by their creation time in seconds since the epoc. This should ensure unique filenames
G.graph_properties['ID'] = G.new_graph_property('int64_t',val=int(time.time()*1000))
G.graph_properties['ntype'] = G.new_graph_property('string',val=network_type)
G.graph_properties['probability'] = G.new_graph_property('double',probabilities[0])

G.save(f'{network_path}/{G.gp.ID}.gt')

print(f'analysis of {G.gp.ID}')
## Run the network analyzer function on all the networks
G = get_local_clutsering(G)
G = get_transitivity(G)
G = get_kirchhoff_index(G)
G = get_ave_shortest_path(G)


## run the liner thershold model on the model, for an appropriate number of times
print(f'Running LTM on {G.gp.ID}')
## Check if any seeds existes, and compute the LTM such that a sufficient number of unique_seeds was used

seeds = []
all_spreads = []
new_picks = 1
new_seeds = random.sample(list(np.delete(G.get_vertices(),seeds)),new_picks)

#Running of the LTM on the newly selected seeds
print('#### new runs p={},k={} ###'.format(np.round(G.gp.probability,4),np.max(G.get_out_degrees(new_seeds))))
spreads,_,thresholds_map = linear_threshold_memory_model(G,thresholds,seed_nodes=[new_seeds],max_iter=nodes[0],init_spread=True)

all_spreads += gt.ungroup_vector_property(spreads,range(len(thresholds)))
seeds.append(new_seeds)
# Merge the spreading vectors
new_spreads = gt.group_vector_property(all_spreads)
# Add new properties ot the graph and save
G.vp['infected_step'] = new_spreads
G.gp['thresholds'] = thresholds_map
G.gp['seeds'] = G.new_gp(value_type='python::object',val=seeds)

#### Calculate polarization speeds and save them to seed vertices ###
## fetch if exists; othterwise create
print('calculating polarization speeds')
if 'polarization_speed' not in G.vp:
    polarization_speeds = np.full((G.num_vertices(),len(cascades)*len(thresholds)),-1.0)
    G.vp['polarization_speed'] = G.new_vp(value_type='vector<double>',vals=polarization_speeds)

polarization_speeds = G.vp.polarization_speed.get_2d_array(range(len(cascades)*len(thresholds))).T

for idy,s in enumerate(seeds):
    # Hack to get around index [0] or [0][0]
    first = polarization_speeds[s]
    while isinstance(first,np.ndarray):
        first = first[0]
    if  first != -20:
        # print('calc')
        #initialize lits to aggregate the polarization speeds
        speeds = []
        # load simulation data
        spread = G.vp.infected_step.get_2d_array(np.arange(idy*len(thresholds),(idy+1)*len(thresholds)))
        for idx,th in enumerate(G.gp.thresholds):
            cascade_sizes = cascades ## Copy where we can shorten list on iterations in while loop
            infected = 0             ## Reset on loop start

            ## Find nodes infected at given step
            val,counts = np.unique(spread[idx],return_counts=True)
            ## Refactor to fraction of nodes
            counts = counts / G.num_vertices()
            ## For each step in LTM add newly infected to total
            for i,new in enumerate(counts[val>-2]):
                infected += new
                ## Current number of infected are used to find polarization speed if larger than cascade size, and the step is after initialization, and all cascade sizes has not been exceeded yet.
                while len(cascade_sizes) > 0 and infected > cascade_sizes[0] and val[i] > 0:
                    speeds.append(infected/val[i])
                    cascade_sizes = cascade_sizes[1:]
            for j in cascade_sizes:
                speeds.append(-1)
        polarization_speeds[s]=speeds

G.vp['polarization_speed'] = G.new_vp(value_type='vector<double>',vals=polarization_speeds)
G.gp['cascades'] = G.new_gp(value_type='python::object',val=cascades)
G.save(f'{network_path}/{G.gp.ID}.gt')
