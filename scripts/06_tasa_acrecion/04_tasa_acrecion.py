"""
Calcula la tasa de acrecion de masa (Mdot_acc) a partir de W10% (Halpha),
usando la calibracion de Natta, Testi & Muzerolle (2004), y genera las
figuras de comparacion con masa, temperatura efectiva y edad.

Requiere: fulldisk.head, RESULTADOS_W10,
          KINCANDS_2MASS_W4_or_M24_fulldisk.cat,
          catalogo_masas_edades_149.csv (script 03)

Genera:
- w10_por_estrella.csv
- catalogo_acrecion.csv
- fig_hist_macc.png, fig_macc_masa.png, fig_macc_teff.png, fig_macc_edad.png
"""
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = './'

# -----------------------------------------------------------
# 1. Correlacionar cada espectro (fulldisk.head) con su
#    gaia_source_id mediante posicion celeste
# -----------------------------------------------------------
head = pd.read_csv(DATA_DIR + 'fulldisk.head', sep=r'\s+', header=None,
                    names=['archivo', 'obsid', 'RA_sex', 'DEC_sex'])

coords = SkyCoord(ra=head['RA_sex'], dec=head['DEC_sex'], unit=(u.hourangle, u.deg))
head['RA_deg'] = coords.ra.deg
head['DEC_deg'] = coords.dec.deg

with open(DATA_DIR + 'KINCANDS_2MASS_W4_or_M24_fulldisk.cat') as f:
    header_line = f.readline()
bcols = header_line.lstrip('#').split()
base = pd.read_csv(DATA_DIR + 'KINCANDS_2MASS_W4_or_M24_fulldisk.cat', sep=r'\s+',
                    skiprows=1, header=None, names=bcols, quotechar='"')

base_coords = SkyCoord(ra=base['ra_epoch2000'].values * u.deg, dec=base['dec_epoch2000'].values * u.deg)
head_coords = SkyCoord(ra=head['RA_deg'].values * u.deg, dec=head['DEC_deg'].values * u.deg)

idx, sep2d, _ = head_coords.match_to_catalog_sky(base_coords)
head['gaia_source_id'] = base['source_id_1'].values[idx]
head['sep_arcsec'] = sep2d.arcsec
print('Separacion maxima (arcsec):', head['sep_arcsec'].max())

# -----------------------------------------------------------
# 2. Combinar W10% por estrella y calcular Mdot_acc
# -----------------------------------------------------------
w10 = pd.read_csv(DATA_DIR + 'RESULTADOS_W10')
merged = w10.merge(head[['archivo', 'gaia_source_id']], on='archivo', how='left')

grouped = merged.groupby('gaia_source_id').agg(
    W10_mean=('W10_km_s', 'mean'),
    W10_std=('W10_km_s', 'std'),
    n_spec=('W10_km_s', 'count')
).reset_index()
grouped['W10_std'] = grouped['W10_std'].fillna(0)
grouped.to_csv('w10_por_estrella.csv', index=False)

# Ecuacion de Natta, Testi & Muzerolle (2004)
grouped['log_Mdot_acc'] = -12.89 + 9.7e-3 * grouped['W10_mean']

# Error riguroso: incertidumbre de los coeficientes + dispersion de W10%
# entre espectros de una misma estrella, sumadas en cuadratura
grouped['e_log_Mdot_acc'] = np.sqrt(
    0.30 ** 2 + (0.7e-3 * grouped['W10_mean']) ** 2 + (9.7e-3 * grouped['W10_std']) ** 2
)
grouped['acretora'] = grouped['W10_mean'] > 200  # umbral de calibracion Natta+04

acretoras = grouped[grouped['acretora']].copy()
print('Estrellas con W10 > 200 km/s:', len(acretoras))

# -----------------------------------------------------------
# 3. Cruzar con el catalogo de masa/Teff/edad de MassAge
# -----------------------------------------------------------
props = pd.read_csv(DATA_DIR + 'catalogo_masas_edades_149.csv')
df = acretoras.merge(props, left_on='gaia_source_id', right_on='Source', how='inner')
df.to_csv('catalogo_acrecion.csv', index=False)
print('Con masa/Teff disponibles:', len(df))
print(df['log_Mdot_acc'].describe())

# -----------------------------------------------------------
# 4. Figuras
# -----------------------------------------------------------
good = df[df['qflag'] == 0]
bad = df[df['qflag'] != 0]

# Histograma
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(df['log_Mdot_acc'], bins=15, color='#c0704a', edgecolor='black', alpha=0.85)
ax.axvline(df['log_Mdot_acc'].mean(), color='#2c3e50', linestyle='--', lw=1.5,
           label=f"Media = {df['log_Mdot_acc'].mean():.2f}")
ax.axvline(df['log_Mdot_acc'].median(), color='#27ae60', linestyle='--', lw=1.5,
           label=f"Mediana = {df['log_Mdot_acc'].median():.2f}")
ax.set_xlabel('log(Mdot_acc / M_sun ano^-1)')
ax.set_ylabel('Numero de estrellas')
ax.set_title(f'Distribucion de la tasa de acrecion (N={len(df)})')
ax.legend(fontsize=9); ax.grid(True, ls=':', alpha=0.5)
plt.tight_layout(); plt.savefig('fig_hist_macc.png', dpi=300); plt.close()


def plot_relacion(x_col, xlabel, log_x, invert_x, fname, color, ref_line=False):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(good[x_col], good['log_Mdot_acc'], yerr=good['e_log_Mdot_acc'],
                fmt='o', color=color, ecolor='gray', elinewidth=0.5, markersize=5,
                markeredgecolor='k', markeredgewidth=0.3, alpha=0.85, zorder=3)
    if len(bad):
        ax.scatter(bad[x_col], bad['log_Mdot_acc'], color='red', marker='x', s=100,
                   linewidth=2, label='qflag!=0', zorder=5)
    if ref_line:
        mgrid = np.linspace(df[x_col].min(), df[x_col].max(), 50)
        norm = df['log_Mdot_acc'].median() - 2 * np.log10(df[x_col].median())
        ax.plot(mgrid, norm + 2 * np.log10(mgrid), 'k--', lw=1, alpha=0.6,
                label='Mdot_acc ~ M*^2 (referencia)')
    if log_x: ax.set_xscale('log')
    if invert_x: ax.invert_xaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel('log(Mdot_acc / M_sun ano^-1)')
    if len(bad) or ref_line: ax.legend(fontsize=9)
    ax.grid(True, ls=':', alpha=0.5, which='both')
    plt.tight_layout(); plt.savefig(fname, dpi=300); plt.close()


plot_relacion('mass_median', 'Masa (M_sun, escala log)', True, False,
              'fig_macc_masa.png', '#4a7c9b', ref_line=True)
plot_relacion('T_median', 'Temperatura efectiva (K)', False, True,
              'fig_macc_teff.png', '#9b6a4a')
plot_relacion('age_median', 'Edad (Ma, escala log)', True, False,
              'fig_macc_edad.png', '#5a8a5a')

print('r(Mdot_acc, masa)  =', np.corrcoef(np.log10(df['mass_median']), df['log_Mdot_acc'])[0, 1])
print('r(Mdot_acc, Teff)  =', np.corrcoef(df['T_median'], df['log_Mdot_acc'])[0, 1])
print('r(Mdot_acc, edad)  =', np.corrcoef(np.log10(df['age_median']), df['log_Mdot_acc'])[0, 1])
