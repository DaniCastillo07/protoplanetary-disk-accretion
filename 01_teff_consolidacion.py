"""
Consolidacion de la temperatura efectiva (Teff) para las 149 estrellas de la
muestra, a partir de las 4 fuentes jerarquicas: SpTClass (P&M13), GBM,
APOGEEnet y TIC v8.2.

Genera:
- teff_consolidado.csv
- histograma_teff.png
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# -----------------------------------------------------------
# Rutas de entrada (ajustar a tu maquina)
# -----------------------------------------------------------
DATA_DIR = './'  # carpeta donde tienes los 4 archivos fuente

def load(fn):
    with open(DATA_DIR + fn) as f:
        header_line = f.readline()
    cols = header_line.lstrip('#').split()
    return pd.read_csv(DATA_DIR + fn, sep=r'\s+', skiprows=1, header=None, names=cols)

opt = load('fulldisk_optspec.txt')
gbm = load('fulldisk_GBM.txt')
apo = load('fulldisk_apogee.txt')
tic = load('fulldisk_TICv82.txt')

# Nivel 1: SpTClass, escala Pecaut & Mamajek (2013)
t1 = pd.DataFrame({
    'gaia_source_id': opt['gaia_source_id'],
    'Teff': opt['TeffP'],
    'eTeff': opt['TerrP'],
    'fuente': 'SpTClass (P&M13)'
})

# Nivel 2: GBM
t2 = pd.DataFrame({
    'gaia_source_id': gbm['gaia_source_id'],
    'Teff': gbm['Teff_GBM'],
    'eTeff': gbm['eTeff_GBM'],
    'fuente': 'GBM'
})

# Nivel 3: APOGEEnet (logTeff -> Teff lineal)
teff_apo = 10 ** apo['logTeff']
eteff_apo = teff_apo * np.log(10) * apo['e_logTeff']
t3 = pd.DataFrame({
    'gaia_source_id': apo['gaia_source_id'],
    'Teff': teff_apo,
    'eTeff': eteff_apo,
    'fuente': 'APOGEEnet'
})

# Nivel 4: TIC v8.2
t4 = pd.DataFrame({
    'gaia_source_id': tic['gaia_source_id'],
    'Teff': tic['Teff'],
    'eTeff': tic['s_Teff'],
    'fuente': 'TIC v8.2'
})

final = pd.concat([t1, t2, t3, t4], ignore_index=True)
final = final.sort_values('gaia_source_id').reset_index(drop=True)
final.to_csv('teff_consolidado.csv', index=False)

print('Total estrellas con Teff:', len(final))
print(final['fuente'].value_counts())
print(final['Teff'].describe())

# -----------------------------------------------------------
# Histograma
# -----------------------------------------------------------
teff = final['Teff']
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(teff, bins=18, color='#4a7c9b', edgecolor='black', alpha=0.85)
ax.axvline(teff.mean(), color='#c0392b', linestyle='--', lw=1.5,
           label=f'Media = {teff.mean():.0f} K')
ax.axvline(teff.median(), color='#27ae60', linestyle='--', lw=1.5,
           label=f'Mediana = {teff.median():.0f} K')
ax.set_xlabel('Temperatura efectiva (K)')
ax.set_ylabel('Numero de estrellas')
ax.set_title('Distribucion de las temperaturas efectivas de la muestra')
ax.legend(fontsize=9)
ax.grid(True, ls=':', lw=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig('histograma_teff.png', dpi=300)
