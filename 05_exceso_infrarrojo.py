"""
Relaciona la tasa de acrecion (Mdot_acc) con el exceso infrarrojo,
cuantificado mediante los indices K-W3 y K-[24] ya definidos en la
metodologia de seleccion de la muestra (Seccion 4.4-4.6).

Requiere: KINCANDS_2MASS_W4_or_M24_fulldisk.cat, catalogo_acrecion.csv (script 04)

Genera:
- catalogo_acrecion_IR.csv
- fig_macc_ir_excess.png
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = './'

with open(DATA_DIR + 'KINCANDS_2MASS_W4_or_M24_fulldisk.cat') as f:
    header_line = f.readline()
bcols = header_line.lstrip('#').split()
base = pd.read_csv(DATA_DIR + 'KINCANDS_2MASS_W4_or_M24_fulldisk.cat', sep=r'\s+',
                    skiprows=1, header=None, names=bcols, quotechar='"')

acc = pd.read_csv(DATA_DIR + 'catalogo_acrecion.csv')

df = acc.merge(base[['source_id_1', 'Kmag', 'W3mag', 'Mag24']],
                left_on='gaia_source_id', right_on='source_id_1', how='left')

# Mag24 = 999.9 es el valor centinela de "sin deteccion/cobertura"
df.loc[df['Mag24'] > 50, 'Mag24'] = np.nan

df['K_W3'] = df['Kmag'] - df['W3mag']
df['K_24'] = df['Kmag'] - df['Mag24']
df.to_csv('catalogo_acrecion_IR.csv', index=False)

print('Con K-W3 valido:', df['K_W3'].notna().sum())
print('Con K-24 valido:', df['K_24'].notna().sum())

r_w3 = np.corrcoef(df['K_W3'], df['log_Mdot_acc'])[0, 1]
sub24 = df.dropna(subset=['K_24'])
r_24 = np.corrcoef(sub24['K_24'], sub24['log_Mdot_acc'])[0, 1]
print('r (K-W3):', r_w3, ' r (K-24):', r_24)

fig, axes = plt.subplots(1, 2, figsize=(13, 6))

good = df[df['qflag'] == 0]
ax = axes[0]
g = good.dropna(subset=['K_W3'])
ax.errorbar(g['K_W3'], g['log_Mdot_acc'], yerr=g['e_log_Mdot_acc'], fmt='o', color='#8a5a9b',
            ecolor='gray', elinewidth=0.5, markersize=5, markeredgecolor='k',
            markeredgewidth=0.3, alpha=0.85)
ax.axvline(0.50, color='gray', ls=':', lw=1, alpha=0.7)
ax.set_xlabel('K - W3 (mag)'); ax.set_ylabel('log(Mdot_acc / M_sun ano^-1)')
ax.set_title(f'K-W3 (N={len(df)}, r={r_w3:.2f})')
ax.grid(True, ls=':', alpha=0.5)

ax = axes[1]
g24 = good.dropna(subset=['K_24'])
ax.errorbar(g24['K_24'], g24['log_Mdot_acc'], yerr=g24['e_log_Mdot_acc'], fmt='o', color='#5a8a8a',
            ecolor='gray', elinewidth=0.5, markersize=5, markeredgecolor='k',
            markeredgewidth=0.3, alpha=0.85)
ax.set_xlabel('K - [24] (mag)'); ax.set_ylabel('log(Mdot_acc / M_sun ano^-1)')
ax.set_title(f'K-[24] (N={len(sub24)}, r={r_24:.2f})')
ax.grid(True, ls=':', alpha=0.5)

plt.tight_layout()
plt.savefig('fig_macc_ir_excess.png', dpi=300)
