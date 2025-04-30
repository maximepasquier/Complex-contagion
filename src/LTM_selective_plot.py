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
    
Output :
    - plots...

Le script analyse les résultats des simulations de LTM pour l'ensemble des paramètres spécifés et produit des graphiques 
de la vitesse de polarisation en fonction du seuil. Les courbes se regroupent suivant :
    - Les mêmes réseaux
    - Les mêmes cascades
    - Des paramètres de mémoire différents (Tau et Persuasion)
'''  

#* Paramètres de simulation
cascade = 0.3 # select only 30% cascades

#* Curves to plot
curve_list = []
#curve_i = (<network_root>,<network_class>,<N>,<K>,<Tau>,<Persuasion>)
'''
curve_list.append(("LTM","ws",1000,8,2,0))
curve_list.append(("LTM","ws",1000,8,4,0))
curve_list.append(("LTM","ws",1000,8,16,0))
curve_list.append(("LTM","ws",1000,8,64,0))

curve_list.append(("LTM","ws",1000,16,2,0))
curve_list.append(("LTM","ws",1000,16,4,0))
curve_list.append(("LTM","ws",1000,16,16,0))
curve_list.append(("LTM","ws",1000,16,64,0))

curve_list.append(("LTM","ws",1000,8,0,0.01))
curve_list.append(("LTM","ws",1000,8,0,0.1))
curve_list.append(("LTM","ws",1000,8,0,0.5))
curve_list.append(("LTM","ws",1000,8,0,1))

curve_list.append(("LTM","ws",1000,16,0,0.01))
curve_list.append(("LTM","ws",1000,16,0,0.1))
curve_list.append(("LTM","ws",1000,16,0,0.5))
curve_list.append(("LTM","ws",1000,16,0,1))

curve_list.append(("LTM","ws",1000,16,32,0.1))
curve_list.append(("LTM","ws",1000,16,32,0.5))
curve_list.append(("LTM","ws",1000,16,64,0.1))
curve_list.append(("LTM","ws",1000,16,64,0.5))
'''
curve_list.append(("LTM","ws",1000,16,32,0))
curve_list.append(("LTM_latence","ws",1000,16,32,0))
curve_list.append(("LTM","ws",1000,16,64,0))
curve_list.append(("LTM_latence","ws",1000,16,64,0))

#* Read CSV files and get data
curve_data = []
for curve in curve_list:
    network_root = curve[0]
    network_type = curve[1]
    n_nodes = curve[2]
    neighbor_k = curve[3]
    tau = curve[4]
    pers = curve[5]
    #* Charger les fichiers polarization.csv et props.csv pour la curve
    nets_prop_file =  f'{network_root}/{network_type}/{n_nodes}/{neighbor_k}/{tau}/{pers}/props.csv'
    polarization_file = f'{network_root}/{network_type}/{n_nodes}/{neighbor_k}/{tau}/{pers}/polarization.csv'

    network_props=pd.read_csv(nets_prop_file,sep='\t')
    network_props.set_index(['ID','network','p'],inplace=True)

    polarization=pd.read_csv(polarization_file,sep='\t')
    polarization.set_index(['ID','network','p','th','seed'],inplace=True) # crée 4 niveaux de multiindexes
    
    # Get la liste des probabilités
    probabilities = np.sort(polarization.index.get_level_values(2).unique())

    # joindre en faisant la moyenne de toutes les simulations avec les paramètres p, th et network identiques
    # Ceci fusionne toutes les occurences ou la seed change mais pas les autres paramètres
    # mpol est une dataframe avec comme colonnes les valeurs moyennées pour les niveau de cascades
    # 'p' 'th' et 'network' sont des multiindexes (pas apparent dans la df)
    mpol = polarization.groupby(['p','th','network']).mean()
    #Get mean polarization on network from seeds
    # Similaire à mpol mais avec 'ID' en plus
    polarization_mean = polarization.groupby(level=[0,1,2,3]).mean()
    curve_data.append((mpol,network_props,probabilities,curve))
    
#* Plotting
for idx,p in enumerate(curve_data[0][2][::-1]):
    
    ix = pd.IndexSlice
    
    fig,axs = plt.subplots(figsize=(7,3),sharex=True,sharey=False,tight_layout=True)
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
    axs.get_xaxis().get_major_formatter().set_scientific(False)
    ### Plotting loop
    C_label = str(network_props.loc[ix[:,network_type,p],:]['CC'].mean().round(2))
    T_label = str(network_props.loc[ix[:,network_type,p],:]['T'].mean().round(2))
    l_label = str(network_props.loc[ix[:,network_type,p],:]['SP'].mean().round(1))
    r_label = str((network_props.loc[ix[:,network_type,p],:]['Rg'].mean()/1000).round(1))
    r_label = r_label + r'\! \times \! 10^{3}'
    title_label_string =r'${}  |  {}  |  {}  |  {} $'.format(C_label,T_label,l_label,r_label)
    title = "C   |   T   |   L   |   Rg"
    axs.set_title(f'Cascade size {cascade} \n {title} \n {title_label_string}')
    
    colors = ['darkslateblue','darkcyan','coral','blue']
    #Toggle hatch
    hatching = True
    ## Pick which cascade sizes are consider, valid choises are 0.1,0.2,...,0.9
    #cascades = list(map(str,list(np.round(np.linspace(0.1,0.9,9),1))))

    save = True

    hatch_regions = {'simple':{'x':[-10,0.12],'y1':[-10,-10],'y2':[10,10],'hatch':'////','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10},
    'trans':{'x':[polarization.index.get_level_values(3).unique()[3],polarization.index.get_level_values(3).unique()[8]],'y1':[-10,-10],'y2':[10,10],'hatch':'_','facecolor':'w','alpha':0,'edgecolor':'black','linewidth':1.0,'zorder':-10 },
    'complex':{'x':[polarization.index.get_level_values(3).unique()[8],10],'y1':[-10, -10],'y2':[10,10],'hatch':'\\\\\\\\ ','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10}}

    for data in curve_data:
        mpol,network_props,probabilities,curve = data

        print('plotting')
        
        #For all values
        probabilities_index_list = np.arange(len(probabilities))
        pick_props_ws = {'ws':probabilities_index_list}
        #For the 4 values used in the paper
        #pick_props_mhk = {'mhk':[0,3,6,10]}

        flag_l_label=0
        probabilities = np.sort(mpol.loc[ix[:,:,network_type],:].index.get_level_values(0).unique())[pick_props_ws[network_type]]
                
        # Labels
        tau_label = str(curve[4])
        pers_label = str(curve[5])
        
        label_string =r'${}  |  {}$'.format(tau_label,pers_label)
            
        pol_fig_legend_label = label_string
        
        # Plot les courbes pour les simulations avec mémoire
        axs.plot(mpol.loc[ix[p,:,network_type],f'{cascade}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cascade}'],ls='--', label=pol_fig_legend_label)
                            

        # Set scale stuff
        axs.set_yscale('log')
        axs.set_ylabel(r'Polarization Speed $(v)$',labelpad=2.5)
        axs.set_xlabel(r'Threshold $( \theta )$',labelpad=2,math_fontfamily='cm')
        axs.set_ylim([1*10**-4,1.1])
        axs.set_xlim([-0.02,0.56])
        ## Legend and title
        legend0 = axs.legend(title=r' $  Tau \;\, |\;\;\, Pers  $', framealpha=1, facecolor='white',loc=[1.1,0],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
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

            axs.text(0.055,0.5,r'{\fontfamily{phv}\selectfont  \textbf{Simple}}',transform=axs.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
            axs.text(0.69,0.9,r'{\fontfamily{phv}\selectfont   \textbf{Complex}}',transform=axs.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
        axs.grid(False)
        axs.grid(False)
    # plt.show()
    # fig.savefig(f'figures/fig1/{network}/fig1_{cas}.pdf')
    if save:
        fig_path = f'figs/selective_plot/ws/{n_nodes}/{neighbor_k}'
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + f'/LTM_c={cascade}_p={p}.pdf')
        print('fig saved')
    else:
        fig.show()