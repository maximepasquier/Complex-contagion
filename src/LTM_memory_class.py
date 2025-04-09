from utils import *
from network_functions import *

'''
Classe pour les simulations des LTM.

Input : 
    - network_root : nom du dossier contenant les réseaux + résultats des simulations
    - network_class : liste des classes de réseau
    - N : nombre de noeuds
    - K : nombre de voisins (moyen ??)
    - Tau : nombre d'itérations prises en compte dans le mécanisme d'inertie
    - Persuasion : valeurs d'influence pour la mécanisme de persuasion
    - cascades : liste de valeurs de cascades
    - threshold : liste de valeurs de seuil
    - probabilities : vecteur de probabilités de rewiring (pour de multiples réseaux ws)
    
Output :
    - Génère les graphes en format gt et les sauvergarde
    - Génère les fichiers polarization.csv, props.csv et config.ini
    - Génère les figures pour l'ensemble des paramètres
'''


class LTM_memory:
    def __init__(self,network_root,network_class,N,K,Tau,Persuasion,cascades,threshold,probabilities):
        self.network_root = network_root
        self.network_class = network_class
        self.N = N
        self.K = K
        self.Tau = Tau
        self.Persuasion = Persuasion
        self.cascades = cascades
        self.threshold = threshold
        self.probabilities = probabilities
        self.cols=['ID', 'network', 'p','th', 'seed']+ cascades.astype('str').tolist()
        self.b = {'ID':[],'network':[], 'CC':[],'T':[],'p':[],'SP':[]}

    '''
    def generate(self):
        1 : Boucle sur les combinaisons de tous les paramètres
        2 : Créer les graphes manquants
        3 : Exécuter le modèle LTM
        4 : Ecrire et enregistrer les fichiers polarization.csv, props.csv et config.ini
    '''
    def generate(self):
        for network_type in self.network_class:
            for n in self.N:        
                for k in self.K:
                    for t in self.Tau:
                        for pers in self.Persuasion:
                            # Génère un vecteur de poids pour le mécanisme d'inertie
                            weights = gen_weights(t)
                            # The polirization file counter
                            count = 0
                            nets_prop_file =  f'{self.network_root}/{network_type}/{n}/{k}/{t}/{pers}/props.csv'
                            polarization_file = f'{self.network_root}/{network_type}/{n}/{k}/{t}/{pers}/polarization.csv'
                            ## If the network_path does not exist it is created
                            if not os.path.exists(f'{self.network_root}/{network_type}/{n}/{k}/{t}/{pers}'):
                                os.makedirs(f'{self.network_root}/{network_type}/{n}/{k}/{t}/{pers}')
                                
                            if not os.path.exists(f'{self.network_root}/{network_type}/Networks/n={n}/k={k}'):
                                os.makedirs(f'{self.network_root}/{network_type}/Networks/n={n}/k={k}')
                                
                            config_path = f'{self.network_root}/{network_type}/{n}/{k}/{t}/{pers}/config.ini'

                            config = configparser.ConfigParser()

                            config['Configuration'] = {'cascade': self.cascades, 'threshold' : self.threshold, 'probabilities' : self.probabilities, 'Weights' : weights, 'Nodes' : n, 'Neighbors' : k, 'Memory' : t, 'Network type' : network_type, 'Persuasion' : pers}

                            # Ecrire la configuration dans un fichier
                            with open(config_path, 'w') as configfile:
                                config.write(configfile)

                            network_props = pd.DataFrame(data=self.b)
                            network_props.set_index(['ID','p'],inplace=True)
                            df = pd.DataFrame(columns = self.cols)
                            
                            # Make a dict with the probabilities as keys, to count the available networks
                            graph_realized=dict.fromkeys(self.probabilities,0)
                            
                            # Count the graphs in network path and create missing networks according to  parameters (N,k,p)
                            for graph_path in get_recursive_graph_paths(f'{self.network_root}/{network_type}/Networks/n={n}/k={k}'):
                                g_old = gt.load_graph(str(graph_path))
                                graph_realized[g_old.gp.probability] = 1
                             
                            
                            #### Create missing networks ####
                            for p in self.probabilities:
                                if graph_realized[p] == 1: # déjà un graph créé pour ce rewiring
                                    continue
                                G = ws_network(n,k,p,seed=None)
                                # Add relevant creation properties with the graph.
                                # Name graphs by their creation time in seconds since the epoc. This should ensure unique filenames
                                G.graph_properties['ID'] = G.new_graph_property('int64_t',val=int(time.time()*1000))
                                G.graph_properties['ntype'] = G.new_graph_property('string',val=network_type)
                                G.graph_properties['probability'] = G.new_graph_property('double',p)
                                G.graph_properties['cascades'] = G.new_gp(value_type='vector<double>',val=self.cascades)

                                G.save(f'{self.network_root}/{network_type}/Networks/n={n}/k={k}/{G.gp.ID}.gt')
                                
                            #### Charger les graphs et lancer les simulations ####
                            for graph_path in list(get_recursive_graph_paths(f'{self.network_root}/{network_type}/Networks/n={n}/k={k}'))[:]:
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
                                    #* Modèle complet (inertie + persuasion)
                                    infected_vectormap, _, _ = linear_threshold_memory_model(G,self.threshold,pers,t,weights,seed_nodes=[seed])
                                    
                                    spread = gt.ungroup_vector_property(infected_vectormap,range(len(self.threshold)))
                                    data = np.empty((len(self.threshold),len(self.cascades),)) * np.nan
                                    for idx,th in enumerate(self.threshold):
                                        speeds = []
                                        # speeds = np.empty((1,len(cascades),)) * np.nan
                                        cascade_sizes = self.cascades ## Copy where we can shorten list on iterations in while loop
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
                                        nan_padding = len(self.cascades) - len(speeds)
                                        speeds = np.pad(speeds, (0, nan_padding), constant_values=np.nan)
                                        # print(speeds)
                                        data[idx,:] = speeds
                                        df.loc[count] = [G.gp.ID] + [G.gp.ntype] + [G.gp.probability] + [th] + [seed] + list(speeds)
                                        count += 1
                                        
                            df.to_csv(polarization_file, sep='\t',index = False)
                            network_props.to_csv(nets_prop_file,sep='\t',mode='w',header=True)
                    
    def visualize(self):
        for network_type in self.network_class:               
            for n_nodes in self.N:
                for neighbor_k in self.K:
                    for tau in self.Tau:
                        for pers in self.Persuasion:
                            ix = pd.IndexSlice
                            colors = ['darkslateblue','darkcyan','coral','blue']
                            #Toggle hatch
                            hatching = True
                            ## Pick which cascade sizes are consider, valid choises are 0.1,0.2,...,0.9
                            cascades = list(map(str,list(np.round(np.linspace(0.1,0.9,9),1))))

                            nets_prop_file =  f'{self.network_root}/{network_type}/{n_nodes}/{neighbor_k}/{tau}/{pers}/props.csv'
                            polarization_file = f'{self.network_root}/{network_type}/{n_nodes}/{neighbor_k}/{tau}/{pers}/polarization.csv'


                            network_props=pd.read_csv(nets_prop_file,sep='\t')
                            network_props.set_index(['ID','network','p'],inplace=True)

                            polarization=pd.read_csv(polarization_file,sep='\t')
                            polarization.set_index(['ID','network','p','th','seed'],inplace=True)


                            mpol = polarization.groupby(['p','th','network']).mean()
                            #Get mean polarization on network from seeds
                            polarization_mean = polarization.groupby(level=[0,1,2,3]).mean()

                            print('plotting')
                            save = True

                            hatch_regions = {'simple':{'x':[-10,0.12],'y1':[-10,-10],'y2':[10,10],'hatch':'////','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10},
                            'trans':{'x':[polarization.index.get_level_values(3).unique()[3],polarization.index.get_level_values(3).unique()[8]],'y1':[-10,-10],'y2':[10,10],'hatch':'_','facecolor':'w','alpha':0,'edgecolor':'black','linewidth':1.0,'zorder':-10 },
                            'complex':{'x':[polarization.index.get_level_values(3).unique()[8],10],'y1':[-10, -10],'y2':[10,10],'hatch':'\\\\\\\\ ','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10}}

                            #For all values
                            pick_props_ws = {'ws':[0,1,2,3,4,5,6,7,8,9,10]}
                            #For the 4 values used in the paper
                            pick_props_mhk = {'mhk':[0,3,6,10]}


                            g = {'param1':[],'param2':[],'th':[]}
                            network_corr = pd.DataFrame(data=g)
                            network_corr.set_index(['param1','param2','th'],inplace=True)
                            methods = ['pearson','spearman']
                            flag_l_label=0
                            probabilities = np.sort(mpol.loc[ix[:,:,network_type],:].index.get_level_values(0).unique())[pick_props_ws[network_type]]
                            for cas in cascades[:]:
                                fig,axs = plt.subplots(figsize=(5.5*1.2,2.31*1),ncols=2,nrows=1,sharex=True,sharey=False,tight_layout=True)
                                fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
                                axs[0].get_xaxis().get_major_formatter().set_scientific(False)
                                ### Plotting loop
                                axs[0].set_title(f'Cascade size {cas}')
                                corr_fig_legend_label = [r'$C$',r'$T$',r'$\ell$',r'$R_g$']
                                for idx,p in enumerate(probabilities[::-1]):
                                    C_label = str(network_props.loc[ix[:,network_type,p],:]['CC'].mean().round(2))
                                    if len(C_label) < 4:
                                        C_label =  C_label + '0'
                                        
                                    T_label = str(network_props.loc[ix[:,network_type,p],:]['T'].mean().round(2))
                                    if len(T_label) < 4:
                                        T_label =  T_label + '0'
                                    
                                    l_label = str(network_props.loc[ix[:,network_type,p],:]['SP'].mean().round(1))
                                    if len(l_label) > 3:
                                        l_label = l_label[:3]
                                        flag_l_label=1
                                    # print(len(l_label))

                                    r_label = str((network_props.loc[ix[:,network_type,p],:]['Rg'].mean()/1000).round(1))
                                    if len(r_label) > 3:
                                        r_label = r_label[:2]
                                    r_label = r_label + r'\! \times \! 10^{3}'

                                    if flag_l_label:
                                        label_string =r'${}  |  {}  |  {}  |  {} $'.format(C_label,T_label,l_label,r_label)
                                    else:
                                        label_string =r'${}  |  {}  |  {} \, |  {} $'.format(C_label,T_label,l_label,r_label)


                                    pol_fig_legend_label = label_string
                                    #Make lines distingushable
                                    if idx < 4:
                                        axs[0].plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],ls='-',label=pol_fig_legend_label)
                                    elif idx < 8:
                                        axs[0].plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],ls='-.',label=pol_fig_legend_label)
                                    elif idx < 12:
                                        axs[0].plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],ls='--',label=pol_fig_legend_label)

                                # Set scale stuff
                                axs[0].set_yscale('log')
                                axs[0].set_ylabel(r'Polarization Speed $(v)$',labelpad=2.5)
                                axs[0].set_xlabel(r'Threshold $( \theta )$',labelpad=2,math_fontfamily='cm')
                                axs[0].set_ylim([1*10**-2,1.1])
                                axs[0].set_xlim([-0.02,0.56])
                                ## Legend and title
                                legend0 = axs[0].legend(title=r' $  C \;\, |\;\;\, T \;\,| \;\,\: \ell  \;\:\,  |\;\;\, R_{g} $', framealpha=1, facecolor='white',loc=[1.1,0],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
                                # legend0 = axs[0].legend(title=r' $  C \;\, | \;\, \ell  \;\:  |\;\;\, R_{g} $', framealpha=1, facecolor='white',loc=[0.535,0.3],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
                                legend0.get_title().set_position((1.5,0))
                                legend0.get_title().set_fontsize('7')
                                for th in polarization.index.get_level_values(3).unique():
                                    (network_corr.loc[('CC','p',th),'speraman'],network_corr.loc[('T','p',th),'speraman'],network_corr.loc[('SP','p',th),'speraman'],network_corr.loc[('Rg','p',th),'speraman'])   = network_props.corrwith((polarization_mean.loc[ix[:,network_type,:,th],f'{cas}'].dropna().groupby('ID').mean()),method=methods[1])[['CC','T','SP','Rg']]


                                for idx,param in enumerate(network_corr.index.get_level_values(0).unique()):
                                    axs[1].plot(network_corr.loc[ix[param,:,:],:].index.get_level_values(2),network_corr.loc[ix[param,:,:],:],label=corr_fig_legend_label[idx],marker=['<','^','>','^'][idx%5],zorder=[3,2,1,4,5][idx%5],c=colors[idx%5],linewidth=2,markersize=[5,5,5,5][idx%5],ls=[':','--','-','--'][idx%5],markeredgewidth=[1,1,1,1][idx%5],markerfacecolor='none')

                                ## set scale stuff
                                axs[1].set_ylabel(r'Correlation $r_s(v,\chi)$',labelpad=2.5)
                                axs[1].set_xlabel(r'Threshold $( \theta )$',labelpad=2)
                                axs[1].set_ylim([-1.05,1.475])
                                axs[1].set_xlim([-0.02,0.51])

                                axs[0].text(-0.02,1.03,r'\textbf{(a)}',transform=axs[0].transAxes,fontsize=8)
                                axs[1].text(-0.02,1.03,r'\textbf{(b)}',transform=axs[1].transAxes,fontsize=8)
                                legend1 = axs[1].legend(title=r'$\chi$',loc=[1.1,0],borderpad=0.2,markerscale=1,handlelength=2,handletextpad=0.4,fontsize=7)

                                legend1.get_title().set_position((0,0))
                                legend1.get_title().set_fontsize('7')

                                box_props = dict(alpha=1,facecolor='w',linewidth=0,zorder=100000,boxstyle='round',pad=0.5)
                                ##Set background color to transitions
                                if hatching:
                                    for x in hatch_regions:
                                        for ax in axs:
                                            ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
                                            ax.fill_between(**hatch_regions[x])
                                            ax.tick_params(axis='x',labelsize=7)
                                            ax.tick_params(axis='y',labelsize=7)

                                for ax in axs:
                                    ax.text(0.055,0.5,r'{\fontfamily{phv}\selectfont  \textbf{Simple}}',transform=ax.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
                                    ax.text(0.69,0.9,r'{\fontfamily{phv}\selectfont   \textbf{Complex}}',transform=ax.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
                                axs[0].grid(False)
                                axs[1].grid(False)
                                # plt.show()
                                # fig.savefig(f'figures/fig1/{network}/fig1_{cas}.pdf')
                                if save:
                                    fig_path = f'figs/ws/{n_nodes}/{neighbor_k}/{tau}/{pers}'
                                    if not os.path.exists(fig_path):
                                        os.makedirs(fig_path)
                                    fig.savefig(fig_path + f'/LTM_{cas}.pdf')
                                    print('fig saved')
                                else:
                                    fig.show()
                                
    def run(self):
        self.generate()
        self.visualize()
        