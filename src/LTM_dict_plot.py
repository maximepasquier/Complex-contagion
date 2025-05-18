from utils import *
from itertools import islice
# Étape 1 : Charger les données depuis le fichier JSON
with open(f'LTM/ws/1000/16/100/32/waiting_dct_p=Network_p=0.0.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/32/0/waiting_dct_p=Network_p=0.1.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/32/0/waiting_dct_p=Network_p=0.021544346900318832.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/32/0/waiting_dct_p=Network_p=0.0.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/waiting_dct_p=Network_p=0.21544346900318823.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/waiting_dct_p=Network_p=0.1.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/waiting_dct_p=Network_p=0.021544346900318832.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/waiting_dct_p=Network_p=0.0.json', 'r') as f:
    waiting_dict = json.load(f)
    
#selected_thresholds = list(waiting_dict.keys())[:2] # select keys

#waiting_dict = dict(islice(waiting_dict.items(), 8,12)) # select keys


# Étape 2 : Plot des valeurs pour chaque threshold
plt.figure(figsize=(10, 6))  # Taille de la figure

# Tracer chaque liste associée à un threshold
for threshold, values in waiting_dict.items():
    # Chaque "values" est une liste de nombres
    # Tracer la courbe pour ce threshold
    plt.plot(values, label=f'Threshold {threshold}')

# Étape 3 : Ajouter des labels et une légende
plt.xlabel('Index des valeurs')
plt.ylabel('Valeurs')
plt.title('Tracé des valeurs pour chaque Threshold')
plt.legend()

# Afficher le graphique
plt.show()

# Tracer un bar plot pour chaque threshold
# Paramètres pour le plot
plt.figure(figsize=(25, 6))  # Taille de la figure
bar_width = 0.05  # Largeur des barres
index = range(len(waiting_dict))  # Index des thresholds

# Espacer les barres côte à côte
for i, (threshold, values) in enumerate(waiting_dict.items()):
    # Compter les occurrences de chaque valeur (distribution des valeurs)
    value_counts = {val: values.count(val) for val in set(values)}

    # Extraire les valeurs uniques et leurs comptes
    unique_values = list(value_counts.keys())
    counts = list(value_counts.values())
    
    # Décaler les positions des barres pour chaque threshold
    plt.bar([x + i * bar_width for x in range(len(unique_values))], counts, bar_width, alpha=0.6, label=f'Threshold {threshold}')

# Étape : Ajouter des labels et une légende
plt.xlabel('Valeurs')
plt.ylabel('Fréquence')
plt.title('Distribution des valeurs pour tous les Thresholds')
plt.legend()

# Afficher le graphique
plt.show()