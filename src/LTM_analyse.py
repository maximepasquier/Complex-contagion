from utils import *
from network_functions import *

'''
Input :
    - network_root : nom du dossier contenant les réseaux
    - network_class : liste des classes de réseau
    - N : liste de nombre de noeuds
    - K : liste de nombre de voisins moyen
    - Tau : liste de taille de la mémoire
    - seed_types : liste de types de seed_nodes
    - fraction_of_seeds : fraction de noeuds pour seed_nodes
    - cascades : liste de taille de cascade
    
Output :
    - /figs/analyse/... : dossier contenant les figures d'analyse des simulations LTM

Le script analyse les résultats des simulations de LTM pour l'ensemble des paramètres spécifés et produit des graphiques 
de la vitesse de polarisation en fonction du seuil. Les courbes se regroupent suivant :
    - Les mêmes réseaux
    - Les mêmes cascades
    - Les mêmes types de seed_nodes
    - Des paramètres de mémoire différents (Tau)
'''  

#* Paramètres de simulation
network_root = 'LTM'
network_class = ['ws','mhk']
N = [1000]
K = [16]
#waiting_counter_max = [0] # nombre max d'itérations admise entre deux pas d'évolution
Tau = [0,1,2,4,8,16,32,64]
seed_types = ['top','rand','bot'] # types de seed_nodes
fraction_of_seeds = 0.05 # fraction de noeuds pour seed_nodes
#Persuasion = [0,0.01,0.05,0.1,0.2,0.5,1]
cascades = np.round(np.linspace(0.1,0.9,9),1)
cascades = [cascades[2],cascades[7]] # select only 30% cascades and 80% cascades

#* Plotting
for network_type in network_class:               
    for n_nodes in N:
        for neighbor_k in K:
            for seed_type in seed_types:
                #for w in waiting_counter_max:
                ix = pd.IndexSlice
                colors = ['darkslateblue','darkcyan','coral','blue']
                #Toggle hatch
                hatching = True
                ## Pick which cascade sizes are consider, valid choises are 0.1,0.2,...,0.9
                #cascades = list(map(str,list(np.round(np.linspace(0.1,0.9,9),1))))

                #* Charger les fichiers polarization.csv et props.csv pour les simulations n'utilisant pas les mécanismes de mémoire (props.csv reste le même dans les 2 cas)
                nets_prop_file =  f'{network_root}/{network_type}/{n_nodes}/{neighbor_k}/{0}/props.csv'
                base_polarization_file = f'{network_root}/{network_type}/{n_nodes}/{neighbor_k}/{0}/{seed_type}{fraction_of_seeds}polarization.csv'

                network_props=pd.read_csv(nets_prop_file,sep='\t')
                network_props.set_index(['ID','network','p'],inplace=True)

                base_polarization=pd.read_csv(base_polarization_file,sep='\t')
                base_polarization.set_index(['ID','network','p','th','seed'],inplace=True) # crée 4 niveaux de multiindexes
                
                # Get la liste des probabilités
                probabilities = np.sort(base_polarization.index.get_level_values(2).unique())

                # joindre en faisant la moyenne de toutes les simulations avec les paramètres p, th et network identiques
                # Ceci fusionne toutes les occurences ou la seed change mais pas les autres paramètres
                # mpol est une dataframe avec comme colonnes les valeurs moyennées pour les niveau de cascades
                # 'p' 'th' et 'network' sont des multiindexes (pas apparent dans la df)
                base_mpol = base_polarization.groupby(['p','th','network']).mean()
                # Get mean polarization on network from seeds
                # Similaire à mpol mais avec 'ID' en plus
                polarization_mean = base_polarization.groupby(level=[0,1,2,3]).mean()

                print('plotting')
                save = True

                hatch_regions = {'simple':{'x':[-10,0.12],'y1':[-10,-10],'y2':[10,10],'hatch':'////','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10},
                'trans':{'x':[base_polarization.index.get_level_values(3).unique()[3],base_polarization.index.get_level_values(3).unique()[8]],'y1':[-10,-10],'y2':[10,10],'hatch':'_','facecolor':'w','alpha':0,'edgecolor':'black','linewidth':1.0,'zorder':-10 },
                'complex':{'x':[base_polarization.index.get_level_values(3).unique()[8],10],'y1':[-10, -10],'y2':[10,10],'hatch':'\\\\\\\\ ','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10}}

                #For all values
                probabilities_index_list = np.arange(len(probabilities))
                if network_type == 'ws':
                    pick_props = {'ws':probabilities_index_list}
                #For the 4 values used in the paper
                elif network_type == 'mhk':
                    pick_props = {'mhk':probabilities_index_list}

                flag_l_label=0
                probabilities = np.sort(base_mpol.loc[ix[:,:,network_type],:].index.get_level_values(0).unique())[pick_props[network_type]]
                        
                for cas in cascades[:]:
                    for idx,p in enumerate(probabilities[::-1]):
                        fig,axs = plt.subplots(figsize=(5,4),sharex=True,sharey=False,tight_layout=True)
                        fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
                        axs.get_xaxis().get_major_formatter().set_scientific(False)
                        ### Plotting loop
                        C_label = str(network_props.loc[ix[:,network_type,p],:]['CC'].mean().round(2))
                        T_label = str(network_props.loc[ix[:,network_type,p],:]['T'].mean().round(2))
                        l_label = str(network_props.loc[ix[:,network_type,p],:]['SP'].mean().round(1))
                        r_label = str((network_props.loc[ix[:,network_type,p],:]['Rg'].mean()/1000).round(1))
                        r_label = r_label + r'\! \times \! 10^{3}'
                        axs.set_title(f'Cascade size {int(100*cas)}\%, {seed_type} {int(100*fraction_of_seeds)}\%\n' + r'$Network = {}, N = {}, K = {}$'.format(network_type,n_nodes,neighbor_k))
                        #corr_fig_legend_label = [r'$C$',r'$T$',r'$\ell$',r'$R_g$']
                        
                        for tau in Tau:
                            #for pers in Persuasion:
                            if tau == 0:
                                # Plot la courbe pour la simulation sans mémoire
                                axs.plot(base_mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),base_mpol.loc[ix[p,:,network_type],f'{cas}'],ls='-',label=r'$0$')
                                continue
                            #* Charger le fichier polarization.csv pour les simulations utilisant les mécanismes de mémoire
                            eval_polarization_file = f'{network_root}/{network_type}/{n_nodes}/{neighbor_k}/{tau}/{seed_type}{fraction_of_seeds}polarization.csv'
                            
                            eval_polarization=pd.read_csv(eval_polarization_file,sep='\t')
                            eval_polarization.set_index(['ID','network','p','th','seed'],inplace=True) # crée 4 niveaux de multiindexes
                            
                            eval_mpol = eval_polarization.groupby(['p','th','network']).mean()
                            
                            # Labels
                            tau_label = str(tau)
                            #pers_label = str(pers)
                            
                            label_string =r'${}$'.format(tau_label)
                                
                            pol_fig_legend_label = label_string
                            
                            # Plot les courbes pour les simulations avec mémoire
                            axs.plot(eval_mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),eval_mpol.loc[ix[p,:,network_type],f'{cas}'],ls='-', label=pol_fig_legend_label)
                            
                        # Set scale stuff
                        axs.set_yscale('log')
                        axs.set_ylabel(r'Polarization Speed $(v)$',labelpad=2.5)
                        axs.set_xlabel(r'Threshold $( \theta )$',labelpad=2,math_fontfamily='cm')
                        axs.set_ylim([1*10**-4,1.1])
                        axs.set_xlim([-0.02,0.56])
                        ## Legend and title
                        legend0 = axs.legend(title=r' $  Tau  $', framealpha=1, facecolor='white',loc=[1.1,0],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
                        # legend0 = axs.legend(title=r' $  C \;\, | \;\, \ell  \;\:  |\;\;\, R_{g} $', framealpha=1, facecolor='white',loc=[0.535,0.3],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
                        legend0.get_title().set_position((1.5,0))
                        legend0.get_title().set_fontsize('7')

                        box_props = dict(alpha=1,facecolor='w',linewidth=0,zorder=100000,boxstyle='round',pad=0.5)
                        ##Set background color to transitions
                        if hatching:
                            for x in hatch_regions:
                                    axs.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
                                    axs.fill_between(**hatch_regions[x])
                                    axs.tick_params(axis='x',labelsize=7)
                                    axs.tick_params(axis='y',labelsize=7)

                            axs.text(0.07,0.1,r'{\fontfamily{phv}\selectfont  \textbf{Simple}}',transform=axs.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
                            axs.text(0.69,0.9,r'{\fontfamily{phv}\selectfont   \textbf{Complex}}',transform=axs.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
                        axs.grid(False)
                        axs.grid(False)
                        
                        # Zone de texte pour infos sur le réseau
                        axs.text(0.75, 0.85,"C   |   T   |   L   |   Rg\n" + r'${}  |  {}  |  {}  |  {} $'.format(C_label,T_label,l_label,r_label), transform=axs.transAxes,bbox=box_props,fontsize=8,fontdict={'family':'sans-serif'}, verticalalignment='top', horizontalalignment='center')
                        
                        # plt.show()
                        # fig.savefig(f'figures/fig1/{network}/fig1_{cas}.pdf')
                        if save:
                            fig_path = f'figs/analyse/{network_type}/{n_nodes}/{neighbor_k}'
                            if not os.path.exists(fig_path):
                                os.makedirs(fig_path)
                            fig.savefig(fig_path + f'/Analyse_{network_type}_{n_nodes}_{neighbor_k}_{seed_type}_{fraction_of_seeds}_c={cas}_p={p}.pdf')
                            print('fig saved')
                        else:
                            fig.show()