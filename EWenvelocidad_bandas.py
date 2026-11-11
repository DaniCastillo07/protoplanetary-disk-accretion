# ==============================
# EW Hα por método de bandas
# en λ y en velocidad
# ==============================

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
import pandas as pd
import os

# -----------------------------------------
# 1. Bandas (ajustadas)
# -----------------------------------------
line_band = (6555, 6570)
cont_left = (6490, 6510)
cont_right = (6600, 6620)

c = 3e5
lambda_0 = 6562.8

# -----------------------------------------
# 2. Carpeta
# -----------------------------------------
carpeta_graficas = "GRAFICAS_EW_BANDAS"
os.makedirs(carpeta_graficas, exist_ok=True)

# -----------------------------------------
# 3. Función
# -----------------------------------------
def medir_EW_bandas(filename):
    try:
        with fits.open(filename) as hdul:
            flux = hdul[2].data.field('FLUX')[0]
            wl_vac = hdul[2].data.field('WAVELENGTH')[0]

        wave_air = pyasl.vactoair2(wl_vac)

        # -------------------------
        # CONTINUO
        # -------------------------
        mask_left = (wave_air >= cont_left[0]) & (wave_air <= cont_left[1])
        mask_right = (wave_air >= cont_right[0]) & (wave_air <= cont_right[1])

        if np.sum(mask_left) == 0 or np.sum(mask_right) == 0:
            return None

        cont_left_flux = np.median(flux[mask_left])
        cont_right_flux = np.median(flux[mask_right])

        continuum = (cont_left_flux + cont_right_flux) / 2

        if continuum <= 0:
            return None

        flux_norm = flux / continuum

        # -------------------------
        # REGIÓN DE LÍNEA
        # -------------------------
        mask_line = (wave_air >= line_band[0]) & (wave_air <= line_band[1])

        if np.sum(mask_line) < 5:
            return None

        wl_line = wave_air[mask_line]
        fl_line = flux_norm[mask_line]

        # quitar outliers fuertes
        fl_line = np.clip(fl_line, 0, 5)

        # -------------------------
        # EW EN λ
        # -------------------------
        EW_lambda = np.trapz(1 - fl_line, wl_line)

        # -------------------------
        # VELOCIDAD
        # -------------------------
        vel_line = c * (wl_line - lambda_0) / lambda_0

        # -------------------------
        # EW EN VELOCIDAD
        # -------------------------
        EW_vel = np.trapz(1 - fl_line, vel_line)
        EW_vel_teo = (c / lambda_0) * EW_lambda

        nombre = os.path.basename(filename)

        print(f"{nombre}")
        print(f"   EW_lambda = {EW_lambda:.3f} Å")
        print(f"   EW_vel = {EW_vel:.3f} km/s")

        # =================================================
        # 📊 GRÁFICA COMBINADA
        # =================================================
        plt.figure(figsize=(10,5))

        plt.plot(wave_air, flux_norm, color='black')
        plt.axhline(1, linestyle='--', color='gray')

        plt.axvspan(cont_left[0], cont_left[1], color='cyan', alpha=0.2)
        plt.axvspan(cont_right[0], cont_right[1], color='cyan', alpha=0.2)
        plt.axvspan(line_band[0], line_band[1], color='red', alpha=0.15)

        plt.xlabel("Wavelength (Å)")
        plt.ylabel("Flujo normalizado")

        ax1 = plt.gca()

        def lambda_to_vel(x):
            return c * (x - lambda_0) / lambda_0

        def vel_to_lambda(v):
            return (v * lambda_0 / c) + lambda_0

        secax = ax1.secondary_xaxis('top', functions=(lambda_to_vel, vel_to_lambda))
        secax.set_xlabel("Velocidad (km/s)")

        plt.title(f"{nombre}\nEW = {EW_lambda:.2f} Å")

        plt.tight_layout()

        ruta = os.path.join(
            carpeta_graficas,
            nombre.replace(".fits", "_bandas.jpg")
        )

        plt.savefig(ruta, dpi=150)
        plt.close()

        return EW_lambda, EW_vel, EW_vel_teo

    except Exception as e:
        print(f"{filename} -> ERROR: {e}")
        return None


# -----------------------------------------
# 4. LOOP
# -----------------------------------------
results = []

carpeta = "FULLDISKS"

for file in os.listdir(carpeta):
    if file.endswith(".fits"):

        ruta = os.path.join(carpeta, file)

        res = medir_EW_bandas(ruta)

        if res is not None:
            EW_lambda, EW_vel, EW_vel_teo = res

            results.append({
                "Archivo": file,
                "EW_lambda_A": EW_lambda,
                "EW_vel_kms": EW_vel,
                "EW_vel_teo_kms": EW_vel_teo
            })


# -----------------------------------------
# 5. CSV
# -----------------------------------------
if results:
    df = pd.DataFrame(results)
    df.to_csv("EW_bandas_comparacion.csv", index=False)
    print("\nCSV guardado: EW_bandas_comparacion.csv")
else:
    print("\nNo hubo resultados válidos")

