"""
Consolida los 3 archivos de resultados de MassAge (cat1_results.csv,
cat2_results.csv, cat3_results.csv) en un unico catalogo de masas y edades.

IMPORTANTE: MassAge guarda el gaia_source_id como float64, lo que produce
perdida de precision en los ultimos digitos de un ID de 19 cifras. Este
script recupera el ID exacto cruzando por RA contra los catalogos de
entrada originales (que si tienen el ID completo).

Requiere: cat1_results.csv, cat2_results.csv, cat3_results.csv,
          massage_cat1_SpT.csv, massage_cat2_lgTeff.csv, massage_cat3_Teff.csv

Genera:
- catalogo_masas_edades_149.csv
- fig_hr_diagram.png, fig_teff_masa.png, fig_edad_masa.png
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = './'


def recover_id(res, inp, tag):
    """Recupera el gaia_source_id de precision completa cruzando por RA
    (redondeado a 6 decimales) contra el catalogo de entrada original."""
    res = res.copy()
    res['RA_r'] = res['RA'].round(6)
    inp = inp.copy()
    inp['RA_r'] = inp['RA'].round(6)
    m = res.merge(inp[['Source', 'RA_r']], on='RA_r', how='left')
    n_missing = m['Source'].isna().sum()
    n_dup = m['Source'].duplicated().sum()
    print(tag, 'filas:', len(res), 'sin match:', n_missing, 'IDs duplicados:', n_dup)
    return m


c1 = pd.read_csv(DATA_DIR + 'cat1_results.csv'); c1.columns = [c.strip() for c in c1.columns]
c2 = pd.read_csv(DATA_DIR + 'cat2_results.csv'); c2.columns = [c.strip() for c in c2.columns]
c3 = pd.read_csv(DATA_DIR + 'cat3_results.csv'); c3.columns = [c.strip() for c in c3.columns]

i1 = pd.read_csv(DATA_DIR + 'massage_cat1_SpT.csv')
i2 = pd.read_csv(DATA_DIR + 'massage_cat2_lgTeff.csv')
i3 = pd.read_csv(DATA_DIR + 'massage_cat3_Teff.csv')

m1 = recover_id(c1, i1, 'cat1'); m1['input_cat'] = 'SpT (optspec)'
m2 = recover_id(c2, i2, 'cat2'); m2['input_cat'] = 'lgTeff (apogee)'
m3 = recover_id(c3, i3, 'cat3'); m3['input_cat'] = 'Teff (GBM+TICv82)'

allc = pd.concat([m1, m2, m3], ignore_index=True)

keep = ['Source', 'RA', 'DEC', 'input_cat', 'Av_wmean', 'eAv_wmean',
        'T_median', 'T_low', 'T_upp', 'T_mad',
        'lgL_median', 'lgL_low', 'lgL_upp', 'lgL_mad',
        'mass_median', 'mass_low', 'mass_upp', 'mass_mad',
        'age_median', 'age_low', 'age_upp', 'age_mad',
        'lgg_median', 'lgg_low', 'lgg_upp', 'lgg_mad',
        'qflag', 'perc_num']
final = allc[keep].copy().sort_values('Source').reset_index(drop=True)
final.to_csv('catalogo_masas_edades_149.csv', index=False)

print('\nTotal:', len(final))
print('qflag distribution:\n', final['qflag'].value_counts())
print('\nMasa (M_sun):\n', final['mass_median'].describe())
print('\nEdad (Ma):\n', final['age_median'].describe())


# -----------------------------------------------------------
# Utilidad de barras de error asimetricas (limite_bajo/alto - mediana)
# -----------------------------------------------------------
def asym_err(df, col):
    lo = (df[f'{col}_median'] - df[f'{col}_low']).clip(lower=0)
    hi = (df[f'{col}_upp'] - df[f'{col}_median']).clip(lower=0)
    return np.vstack([lo, hi])


good = final[final['qflag'] == 0]
bad = final[final['qflag'] != 0]

# --- Diagrama H-R ---
fig, ax = plt.subplots(figsize=(7, 6))
sc = ax.scatter(good['T_median'], good['lgL_median'], c=good['mass_median'],
                 cmap='plasma_r', s=32, edgecolor='k', linewidth=0.3, zorder=3)
ax.errorbar(good['T_median'], good['lgL_median'],
            xerr=asym_err(good, 'T'), yerr=asym_err(good, 'lgL'),
            fmt='none', ecolor='gray', elinewidth=0.4, alpha=0.4, zorder=2)
ax.scatter(bad['T_median'], bad['lgL_median'], color='red', marker='x', s=120,
           linewidth=2.2, label=f'qflag!=0 ({len(bad)} estrellas)', zorder=5)
cbar = plt.colorbar(sc, ax=ax); cbar.set_label('Masa (M_sun)')
ax.invert_xaxis()
ax.set_xlabel('Temperatura efectiva (K)'); ax.set_ylabel('log(L/L_sun)')
ax.set_title('Diagrama H-R de la muestra (149 estrellas)')
ax.grid(True, ls=':', alpha=0.5); ax.legend(fontsize=9, loc='upper right')
plt.tight_layout(); plt.savefig('fig_hr_diagram.png', dpi=300); plt.close()

# --- Teff vs Masa ---
fig, ax = plt.subplots(figsize=(7, 6))
ax.errorbar(good['T_median'], good['mass_median'],
            xerr=asym_err(good, 'T'), yerr=asym_err(good, 'mass'),
            fmt='o', color='#4a7c9b', ecolor='gray', elinewidth=0.5,
            markersize=5, markeredgecolor='k', markeredgewidth=0.3, alpha=0.85, zorder=3)
ax.scatter(bad['T_median'], bad['mass_median'], color='red', marker='x', s=120,
           linewidth=2.2, label=f'qflag!=0 ({len(bad)} estrellas)', zorder=5)
ax.invert_xaxis()
ax.set_xlabel('Temperatura efectiva (K)'); ax.set_ylabel('Masa (M_sun)')
ax.set_title('Relacion temperatura efectiva - masa')
ax.grid(True, ls=':', alpha=0.5); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig('fig_teff_masa.png', dpi=300); plt.close()

# --- Edad vs Masa ---
fig, ax = plt.subplots(figsize=(7, 6))
ax.errorbar(good['age_median'], good['mass_median'],
            xerr=asym_err(good, 'age'), yerr=asym_err(good, 'mass'),
            fmt='o', color='#c0704a', ecolor='gray', elinewidth=0.5,
            markersize=5, markeredgecolor='k', markeredgewidth=0.3, alpha=0.85, zorder=3)
ax.scatter(bad['age_median'], bad['mass_median'], color='red', marker='x', s=120,
           linewidth=2.2, label=f'qflag!=0 ({len(bad)} estrellas)', zorder=5)
ax.set_xscale('log')
ax.set_xlabel('Edad (Ma, escala log)'); ax.set_ylabel('Masa (M_sun)')
ax.set_title('Relacion edad - masa')
ax.grid(True, ls=':', alpha=0.5, which='both'); ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig('fig_edad_masa.png', dpi=300); plt.close()
