import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Lectura de los resultados de los 208 espectros
# ---------------------------------------------------------

df = pd.read_csv("Comparacion_EW_Bandas_Gauss(3).csv")

# ---------------------------------------------------------
# Construccion de intervalos comunes para ambos metodos
# ---------------------------------------------------------

valores = np.concatenate([
    df["EW_Bandas"].values,
    df["EW_Gauss"].values
])

bins = np.histogram_bin_edges(
    valores,
    bins=20
)

# ---------------------------------------------------------
# Generacion del histograma
# ---------------------------------------------------------

fig, ax = plt.subplots(figsize=(8.5, 5.8))

ax.hist(
    df["EW_Bandas"],
    bins=bins,
    alpha=0.55,
    edgecolor="black",
    linewidth=0.7,
    label="Metodo de bandas"
)

ax.hist(
    df["EW_Gauss"],
    bins=bins,
    alpha=0.55,
    edgecolor="black",
    linewidth=0.7,
    label="Ajuste gaussiano"
)

# Linea de referencia en EW = 0
ax.axvline(
    0,
    linestyle="--",
    linewidth=1
)

# ---------------------------------------------------------
# Etiquetas y formato de la figura
# ---------------------------------------------------------

ax.set_xlabel(
    r"Ancho equivalente de H$\alpha$ ($\AA$)"
)

ax.set_ylabel(
    "Numero de espectros"
)

ax.set_title(
    r"Distribucion de los anchos equivalentes de H$\alpha$"
)

ax.legend()

ax.grid(
    axis="y",
    linestyle=":",
    alpha=0.35
)

plt.tight_layout()

# ---------------------------------------------------------
# Exportacion de la figura
# ---------------------------------------------------------

plt.savefig(
    "histograma_EW.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()
