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
    - Les mêmes cascades
    - Des paramètres de mémoire différents (Tau et Persuasion)
    - Des réseaux avec un rewiring différent
'''  

#* Paramètres de simulation
network_root = 'LTM'
network_class = ['ws']
N = [1000]
K = [8,16]
Tau = [0,1,2,4,8,16,32,64]
Persuasion = [0,0.01,0.05,0.1,0.2,0.5,1]
cascades = np.round(np.linspace(0.1,0.9,9),1)
cascades = [cascades[2]]

#* Plotting
for network_type in network_class:               
    for n_nodes in N:
        for neighbor_k in K:
            for tau in Tau:
                for pers in Persuasion:
                    ix = pd.IndexSlice
                    colors = ['darkslateblue','darkcyan','coral','blue']
                    #Toggle hatch
                    hatching = True
                    ## Pick which cascade sizes are consider, valid choises are 0.1,0.2,...,0.9
                    #cascades = list(map(str,list(np.round(np.linspace(0.1,0.9,9),1))))

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

                    print('plotting')
                    save = True

                    hatch_regions = {'simple':{'x':[-10,0.12],'y1':[-10,-10],'y2':[10,10],'hatch':'////','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10},
                    'trans':{'x':[polarization.index.get_level_values(3).unique()[3],polarization.index.get_level_values(3).unique()[8]],'y1':[-10,-10],'y2':[10,10],'hatch':'_','facecolor':'w','alpha':0,'edgecolor':'black','linewidth':1.0,'zorder':-10 },
                    'complex':{'x':[polarization.index.get_level_values(3).unique()[8],10],'y1':[-10, -10],'y2':[10,10],'hatch':'\\\\\\\\ ','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10}}

                    #For all values
                    probabilities_index_list = np.arange(len(probabilities))
                    pick_props_ws = {'ws':probabilities_index_list}
                    #For the 4 values used in the paper
                    #pick_props_mhk = {'mhk':[0,3,6,10]}

                    flag_l_label=0
                    probabilities = np.sort(mpol.loc[ix[:,:,network_type],:].index.get_level_values(0).unique())[pick_props_ws[network_type]]
                    
                    for cas in cascades[:]:
                        fig,axs = plt.subplots(figsize=(7,3),sharex=True,sharey=False,tight_layout=True)
                        fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
                        axs.get_xaxis().get_major_formatter().set_scientific(False)
                        ### Plotting loop
                        axs.set_title(f'Cascade size {cas}')
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
                            if idx == 0:
                                axs.plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],marker='s',markersize = 3,ls='-',label=pol_fig_legend_label)
                            elif idx == 1:
                                axs.plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],marker='^',markersize = 3,ls='-',label=pol_fig_legend_label)
                            elif idx == 2:
                                axs.plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],marker='o',markersize = 3,ls='-',label=pol_fig_legend_label)
                            elif idx == 3:
                                axs.plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],marker='v',markersize = 3,ls='-',label=pol_fig_legend_label)
                            elif idx > 3:
                                axs.plot(mpol.loc[ix[p,:,network_type],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network_type],f'{cas}'],ls='--',label=pol_fig_legend_label)

                        # Set scale stuff
                        axs.set_yscale('log')
                        axs.set_ylabel(r'Polarization Speed $(v)$',labelpad=2.5)
                        axs.set_xlabel(r'Threshold $( \theta )$',labelpad=2,math_fontfamily='cm')
                        axs.set_ylim([1*10**-4,1.1])
                        axs.set_xlim([-0.02,0.56])
                        ## Legend and title
                        legend0 = axs.legend(title=r' $  C \;\, |\;\;\, T \;\,| \;\,\: \ell  \;\:\,  |\;\;\, R_{g} $', framealpha=1, facecolor='white',loc=[1.1,0],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
                        # legend0 = axs[0].legend(title=r' $  C \;\, | \;\, \ell  \;\:  |\;\;\, R_{g} $', framealpha=1, facecolor='white',loc=[0.535,0.3],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
                        legend0.get_title().set_position((1.5,0))
                        legend0.get_title().set_fontsize('7')
                        

                        axs.text(-0.02,1.03,r'\textbf{(a)}',transform=axs.transAxes,fontsize=8)

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
                        # plt.show()
                        # fig.savefig(f'figures/fig1/{network}/fig1_{cas}.pdf')
                        if save:
                            fig_path = f'figs/visualize/ws/{n_nodes}/{neighbor_k}/{tau}/{pers}'
                            if not os.path.exists(fig_path):
                                os.makedirs(fig_path)
                            fig.savefig(fig_path + f'/LTM_{cas}.pdf')
                            print('fig saved')
                        else:
                            fig.show()
