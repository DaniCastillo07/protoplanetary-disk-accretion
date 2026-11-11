import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------------------
# Trayectorias de Hayashi
# Basadas en modelos BHAC15 (Baraffe et al. 2015)
# Formato: (edad_Myr, log_Teff, log_L/Lsun)
# -----------------------------------------
tracks = {
    0.1: [(0.5,3.491,-1.50),(1,3.490,-1.70),(3,3.487,-1.95),
          (10,3.485,-2.15),(30,3.481,-2.30),(100,3.478,-2.50),
          (1000,3.470,-2.80)],
    0.3: [(0.5,3.530,-0.80),(1,3.528,-1.00),(3,3.525,-1.30),
          (10,3.522,-1.55),(30,3.518,-1.75),(100,3.530,-1.90),
          (1000,3.545,-2.10)],
    0.5: [(0.5,3.580,-0.30),(1,3.575,-0.55),(3,3.568,-0.85),
          (10,3.562,-1.10),(30,3.565,-1.25),(100,3.590,-1.35),
          (1000,3.620,-1.50)],
    1.0: [(0.5,3.660, 0.45),(1,3.648, 0.15),(3,3.635,-0.15),
          (10,3.628,-0.38),(30,3.640,-0.40),(100,3.680,-0.30),
          (1000,3.760, 0.00)],
    1.5: [(0.5,3.710, 0.90),(1,3.695, 0.65),(3,3.675, 0.30),
          (10,3.668, 0.10),(30,3.685, 0.12),(100,3.730, 0.25),
          (1000,3.800, 0.50)],
    2.0: [(0.5,3.748, 1.25),(1,3.735, 1.00),(3,3.715, 0.65),
          (10,3.705, 0.45),(30,3.718, 0.48),(100,3.760, 0.60),
          (1000,3.830, 0.85)],
}

# Secuencia principal de edad cero (ZAMS)
zams_teff = np.array([3.46,3.50,3.54,3.58,3.63,3.68,3.72,3.76,
                      3.80,3.85,3.90,3.95,4.00,4.10,4.20,4.40,4.60])
zams_lum  = np.array([-2.8,-2.3,-1.8,-1.3,-0.8,-0.35,0.05,0.35,
                       0.65,1.0,1.35,1.70,2.05,2.60,3.10,4.00,4.80])

isochrone_ages = [1, 3, 10, 30]
iso_colors     = ['#e63946', '#f4a261', '#2a9d8f', '#457b9d']
iso_labels     = ['1 Myr', '3 Myr', '10 Myr', '30 Myr']

masses      = [0.1, 0.3, 0.5, 1.0, 1.5, 2.0]
mass_colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(masses)))

fig, ax = plt.subplots(figsize=(8, 7))
fig.patch.set_facecolor('white')
ax.set_facecolor('#fafafa')

ax.plot(zams_teff, zams_lum, 'k--', lw=1.8, alpha=0.6)

for i, mass in enumerate(masses):
    data  = np.array(tracks[mass])
    teffs = data[:, 1]
    lums  = data[:, 2]
    ax.plot(teffs, lums, '-', color=mass_colors[i], lw=1.8, alpha=0.85)
    ax.annotate(fr'${mass}\,M_\odot$',
                xy=(teffs[0], lums[0]),
                xytext=(4, 3), textcoords='offset points',
                fontsize=8.5, color=mass_colors[i], fontweight='bold')

for age, color, label in zip(isochrone_ages, iso_colors, iso_labels):
    iso_points = []
    for mass in masses:
        data   = np.array(tracks[mass])
        ages_m = data[:, 0]
        if age < ages_m[0] or age > ages_m[-1]:
            continue
        iso_points.append((np.interp(age, ages_m, data[:, 1]),
                           np.interp(age, ages_m, data[:, 2])))
    if iso_points:
        iso_arr = np.array(iso_points)
        ax.plot(iso_arr[:, 0], iso_arr[:, 1], 's--',
                color=color, lw=1.4, ms=5, alpha=0.9, label=label)

ax.set_xlim(3.85, 3.46)
ax.set_ylim(-3.2, 1.8)

teff_ticks  = [3.48,3.52,3.56,3.60,3.65,3.70,3.75,3.80,3.85]
ax.set_xticks(teff_ticks)
ax.set_xticklabels([f'{int(10**t):,}' for t in teff_ticks], fontsize=9)
ax.set_yticks([-3,-2,-1,0,1])
ax.set_yticklabels([r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',
                    r'$10^{0}$',r'$10^{1}$'], fontsize=9)

ax.set_xlabel('Temperatura efectiva (K)', fontsize=11)
ax.set_ylabel(r'Luminosidad ($L/L_\odot$)', fontsize=11)
ax.set_title('Diagrama H-R: Pre-secuencia principal\n'
             '(Estrellas T Tauri)', fontsize=12, pad=10)
ax.grid(True, ls=':', lw=0.6, alpha=0.5)

mass_handles = [Line2D([0],[0], color=mass_colors[i], lw=2,
                label=fr'${m}\,M_\odot$')
                for i, m in enumerate(masses)]
iso_handles  = [Line2D([0],[0], color=c, lw=1.4, ls='--',
                marker='s', ms=5, label=l)
                for c, l in zip(iso_colors, iso_labels)]
zams_handle  = [Line2D([0],[0], color='k', lw=1.8, ls='--',
                alpha=0.6, label='ZAMS')]

leg1 = ax.legend(handles=mass_handles,
                 title='Trayectorias\nde Hayashi',
                 loc='upper right', fontsize=8.5,
                 title_fontsize=9, framealpha=0.9)
ax.legend(handles=iso_handles + zams_handle,
          title='Isocronas', loc='lower left',
          fontsize=8.5, title_fontsize=9, framealpha=0.9)
ax.add_artist(leg1)

ax.axvspan(3.50, 3.75, alpha=0.04, color='orange')
ax.text(3.625, 1.6, 'Region T Tauri', ha='center',
        fontsize=9, color='darkorange', style='italic')
ax.text(0.01, 0.01, 'Basado en Baraffe et al. (2015)',
        transform=ax.transAxes, fontsize=7.5,
        color='gray', va='bottom')

plt.tight_layout()
plt.savefig('hr_diagram.png', dpi=180,
            bbox_inches='tight', facecolor='white')
