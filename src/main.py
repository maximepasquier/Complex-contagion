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
import configparser
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
from LTM_memory_class import *

network_root = 'LTM_networks'
# network_type = 'mhk'
network_class = ['ws']
network_type = network_class[0]

N = [100,200]
K = [4,8]
Tau = [0,2]
Persuasion = [0,0.05]

cascades = np.round(np.linspace(0.1,0.9,9),1)
threshold = np.linspace(0.01,0.5,16)
probabilities = np.append([0], np.logspace(-3,-0,10))

## How many realizations to do of each set of parameters
desired_realizations= 1
## How many unique starting points to run the LTM from on a network
unique_network_seeds = 1

cols=['ID', 'network', 'p','th', 'seed']+ cascades.astype('str').tolist()
b = {'ID':[],'network':[], 'CC':[],'T':[],'p':[],'SP':[]}

#* Simulations

ltm = LTM_memory(network_root,network_class,network_type,N,K,Tau,Persuasion,cascades,threshold,probabilities,desired_realizations,unique_network_seeds,cols,b)

#ltm.visualize()
ltm.run()