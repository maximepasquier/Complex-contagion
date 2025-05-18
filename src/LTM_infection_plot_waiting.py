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
with open(f'LTM/ws/1000/8/1/32/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
    data1 = json.load(f)
with open(f'LTM/ws/1000/8/5/32/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
    data5 = json.load(f)
with open(f'LTM/ws/1000/8/10/32/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
    data10 = json.load(f)
with open(f'LTM/ws/1000/8/20/32/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
    data20 = json.load(f)
with open(f'LTM/ws/1000/8/50/32/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
    data50 = json.load(f)
with open(f'LTM/ws/1000/8/100/32/infected_p=Network_p=0.21544346900318823.json', 'r') as f:
    data100 = json.load(f)

# === Plot ===
plt.figure(figsize=(10, 6))

for data in [data1, data5, data10, data20, data50, data100]:
    for th_str, content in data.items():
        val = content["val"]
        mask = np.array(val) >= -1 # écarter les valeurs négatives par défaut
        counts = content["counts"]
        extended_val, extended_counts = extend_points(np.array(val)[mask],np.cumsum(np.array(counts)[mask]))
        #plt.plot(extended_val, extended_counts, label=f"th={th_str}")
        plt.plot(extended_val, extended_counts)

# === Personnalisation ===
plt.xlabel("val")
plt.ylabel("counts")
plt.title("Superposition des courbes val vs counts")
plt.legend(loc="best", fontsize="small")
plt.grid(True)
plt.tight_layout()

# === Affichage ===
plt.show()
