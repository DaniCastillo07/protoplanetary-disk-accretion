# =====================================================
# W10(Hα) ROBUSTO FINAL
# =====================================================

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.ndimage import gaussian_filter1d


# -----------------------------------------------------
# FUNCIÓN PRINCIPAL
# -----------------------------------------------------

def calcular_W10(wl, fl):

    # -------------------------
    # REGIÓN AMPLIA
    # -------------------------
    mask = (wl >= 6500) & (wl <= 6620)
    wl = wl[mask]
    fl = fl[mask]

    if len(wl) < 20:
        raise ValueError("Muy pocos datos")

    # -------------------------
    # CONTINUO ROBUSTO
    # -------------------------
    ind_izq = (wl >= 6500) & (wl <= 6520)
    ind_der = (wl >= 6600) & (wl <= 6620)

    cont_izq = np.median(fl[ind_izq])
    cont_der = np.median(fl[ind_der])

    continuo = (cont_izq + cont_der) / 2

    if continuo <= 0 or not np.isfinite(continuo):
        raise ValueError("Continuo inválido")

    # -------------------------
    # NORMALIZACIÓN
    # -------------------------
    fl_norm = fl / continuo

    # -------------------------
    # SUAVIZADO (clave)
    # -------------------------
    fl_smooth = gaussian_filter1d(fl_norm, sigma=1)

    # -------------------------
    # PICO Y NIVEL 10%
    # -------------------------
    idx_peak = np.argmax(fl_smooth)
    pico = fl_smooth[idx_peak]

    if pico <= 1:
        raise ValueError("No hay emisión clara")

    nivel_10 = 1 + 0.1 * (pico - 1)

    # -------------------------
    # BUSCAR CRUCES
    # -------------------------
    diff = fl_smooth - nivel_10
    sign = np.sign(diff)

    cambios = np.where(np.diff(sign))[0]

    if len(cambios) < 2:
        raise ValueError("No hay suficientes cruces")

    cruces = []

    for i in cambios:
        x1, x2 = wl[i], wl[i+1]
        y1, y2 = diff[i], diff[i+1]

        x_cross = x1 - y1 * (x2 - x1) / (y2 - y1)
        cruces.append(x_cross)

    cruces = np.array(cruces)

    # -------------------------
    # SOLO CRUCES CERCA DEL PICO
    # -------------------------
    centro_wl = wl[idx_peak]

    cruces_izq = cruces[cruces < centro_wl]
    cruces_der = cruces[cruces > centro_wl]

    if len(cruces_izq) == 0 or len(cruces_der) == 0:
        raise ValueError("Cruces insuficientes cerca del pico")

    wl_izq = cruces_izq[-1]
    wl_der = cruces_der[0]

    W10 = wl_der - wl_izq

    # -------------------------
    # CONVERSIÓN A km/s
    # -------------------------
    c = 3e5
    W10_kms = (W10 / 6562.8) * c

    return W10, W10_kms, wl, fl_norm, fl_smooth, nivel_10, wl_izq, wl_der, centro_wl


# -----------------------------------------------------
# MAIN
# -----------------------------------------------------

if len(sys.argv) < 2:
    print("Uso: python3 ancho_10.py lista.lis")
    sys.exit()

lista_archivos = sys.argv[1]

carpeta_graficas = "GRAFICAS_W10"
carpeta_resultados = "RESULTADOS_W10"

os.makedirs(carpeta_graficas, exist_ok=True)
os.makedirs(carpeta_resultados, exist_ok=True)

with open(lista_archivos) as f:
    archivos = f.read().splitlines()

archivos = [
    j if j.startswith("FULLDISKS/")
    else os.path.join("FULLDISKS", j)
    for j in archivos
]

nombres = []
W10s = []
W10s_kms = []

print("\n===== W10 FINAL =====\n")


# -----------------------------------------------------
# LOOP
# -----------------------------------------------------

for archivo in archivos:

    nombre = os.path.basename(archivo)
    print("Procesando:", nombre)

    try:

        with fits.open(archivo) as hdul:
            fl = hdul[2].data.field('FLUX')[0]
            wl_vac = hdul[2].data.field('WAVELENGTH')[0]

        wl = pyasl.vactoair2(wl_vac, mode='ciddor')

        W10, W10_kms, wl_plot, fl_norm, fl_smooth, nivel_10, wl_izq, wl_der, centro = calcular_W10(wl, fl)

        print(f"   W10 = {W10:.2f} Å  ({W10_kms:.0f} km/s)")

    except Exception as e:
        print("   Error:", e)
        W10 = np.nan
        W10_kms = np.nan

    nombres.append(nombre)
    W10s.append(W10)
    W10s_kms.append(W10_kms)


    # -------------------------------------------------
    # GRÁFICA
    # -------------------------------------------------

    try:
        plt.figure(figsize=(10,5))

        plt.plot(wl_plot, fl_norm, color='gray', alpha=0.5, label="Datos")
        plt.plot(wl_plot, fl_smooth, color='black', label="Suavizado")

        plt.axhline(1, linestyle='--', color='gray', label="Continuo")
        plt.axhline(nivel_10, linestyle='--', color='red', label="Nivel 10%")

        plt.axvline(wl_izq, color='blue', label="W10")
        plt.axvline(wl_der, color='blue')

        plt.axvline(centro, linestyle='--', color='green', label="Centro")

        plt.title(f"{nombre}\nW10(Hα) = {W10:.2f} Å ({W10_kms:.0f} km/s)")
        plt.xlabel("Wavelength (Å)")
        plt.ylabel("Flujo normalizado")

        plt.xlim(6530, 6600)
        plt.grid(alpha=0.2)

        plt.legend()
        plt.tight_layout()

        ruta_fig = os.path.join(
            carpeta_graficas,
            nombre.replace(".fits", "_W10.jpg")
        )

        plt.savefig(ruta_fig, dpi=200)
        plt.close()

    except:
        pass


# -----------------------------------------------------
# GUARDAR CSV
# -----------------------------------------------------

df = pd.DataFrame({
    "archivo": nombres,
    "W10_AA": W10s,
    "W10_km_s": W10s_kms
})

ruta_csv = os.path.join(
    carpeta_resultados,
    lista_archivos.replace(".lis", "_W10.csv")
)

df.to_csv(ruta_csv, index=False)

print("\n===== TERMINADO =====")
print("CSV en:", ruta_csv)
print("Gráficas en:", carpeta_graficas)

