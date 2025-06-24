# Contagion complexe

Les contagions complexes se distinguent des contagions simples par un effet de renforcement social nécessitant l’influence de plusieurs sources. Le Linear Threshold Model (LTM) est un cadre largement étudié pour modéliser ce type de dynamique sociale. Toutefois, aucune étude ne s’est, à notre connaissance, intéressée à l’introduction de mécanismes de mémoire dans ce modèle sans en modifier la logique. Ce travail propose une extension du LTM intégrant une forme d’inertie au changement fondée sur la mémoire des états précédents, et en analyse les effets sur les dynamiques de propagation.

Deux types de pondération temporelle (linéaire et exponentielle) sont explorés, afin de moduler l'influence décroissante du passé sur la prise de décision des agents. Les résultats montrent que l’ajout de mémoire ralentit significativement les dynamiques de contagion et peut, dans certains cas, empêcher l’émergence de cascades globales. L’effet de la mémoire est particulièrement marqué dans les réseaux homogènes de type Watts-Strogatz, tandis que les réseaux hétérogènes de type Modified Holme-Kim y sont moins sensibles.

En parallèle, plusieurs tentatives d’optimisation algorithmique du LTM ont été menées afin de réduire le coût computationnel des simulations sur de grands réseaux. Bien que ces améliorations n’aient pas permis d’obtenir des gains significatifs dans l’environnement Python utilisé, elles ouvrent des pistes prometteuses, notamment pour des implémentations dans des langages plus performants.

Ce travail met en lumière l’intérêt d’introduire des mécanismes de mémoire dans l’étude des contagions sociales et souligne les limites computationnelles actuelles pour simuler efficacement des systèmes de grande taille.

## Utilisation du code

1) Générer les réseaux : LTM_net_gen.py
2) Simuler les LTMs : LTM_simulate.py
3) Produire les résultats : LTM_analyse.py, LTM_selective_analyse.py, LTM_visualization.py, LTM_selective_visualization.py
