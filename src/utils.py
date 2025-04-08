'''
Ensemble des import pour l'exécution de tous les scripts.
Ce fichier est à importer dans chacun des autres scripts.
'''

import numpy as np
import matplotlib.pyplot as plt
import warnings
import networkx as nx
import scipy as sp
from pathlib import Path
import graph_tool as gt
import graph_tool.clustering
import graph_tool.spectral
import graph_tool.topology
from scipy import linalg
import random
from scipy import sparse
import matplotlib.pyplot as plt
import glob
import pandas as pd
import configparser
import os
from scipy.sparse import coo_array
import time
from matplotlib.colors import to_hex
import sys
warnings.filterwarnings("ignore")
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble = r'\usepackage{mathptmx}')