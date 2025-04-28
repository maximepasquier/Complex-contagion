import plotly.graph_objects as go
import numpy as np

x = np.linspace(0, 10, 100)
fig = go.Figure()

# Simulons 10 courbes
for i in range(10):
    y = np.sin(x + i)
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name=f"courbe {i}"))

fig.update_layout(title="Cliquer sur les légendes pour afficher/masquer")
fig.show()
