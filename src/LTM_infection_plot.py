from utils import *

def extend_points(val,counts):
    """
    Fonction pour étendre les points de val et counts
    """
    extended_val = []
    extended_counts = []
    
    tmp = 0
    for v, c in zip(val, counts):
        extended_val.append(v)
        extended_val.append(v)
        extended_counts.append(tmp)
        extended_counts.append(c)
        tmp = c
    
    return extended_val, extended_counts

# === Chargement des données ===
#with open(f'LTM_latence/ws/1000/16/32/0/infected_p=Network_p=0.0.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/32/0/infected_p=Network_p=0.1.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/32/0/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/32/0/infected_p=Network_p=0.021544346900318832.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/infected_p=Network_p=0.0.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/infected_p=Network_p=0.1.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
#with open(f'LTM_latence/ws/1000/16/64/0/infected_p=Network_p=0.021544346900318832.json', 'r') as f:

with open(f'LTM/ws/1000/16/64/0/infected_p=Network_p=0.021544346900318832.json', 'r') as f:
    data = json.load(f)

# === Plot ===
plt.figure(figsize=(10, 6))

for th_str, content in data.items():
    val = content["val"]
    mask = np.array(val) >= -1 # écarter les valeurs négatives par défaut
    counts = content["counts"]
    extended_val, extended_counts = extend_points(np.array(val)[mask],np.cumsum(np.array(counts)[mask]))
    plt.plot(extended_val, extended_counts, label=f"th={th_str}")

# === Personnalisation ===
plt.xlabel("val")
plt.ylabel("counts")
plt.title("Superposition des courbes val vs counts")
plt.legend(loc="best", fontsize="small")
plt.grid(True)
plt.tight_layout()

# === Affichage ===
plt.show()
