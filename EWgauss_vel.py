# =====================================================
# EW(Hα) por ajuste gaussiano (λ + velocidad)
# CON SIGNO FÍSICO CORRECTO (emisión negativa)
# =====================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.optimize import curve_fit

# -----------------------------------------
# CONSTANTES
# -----------------------------------------
c = 3e5
lambda_0 = 6562.8

# -----------------------------------------
# CARPETA
# -----------------------------------------
carpeta_graficas = "GRAFICAS_EW_GAUSS_VEL"
os.makedirs(carpeta_graficas, exist_ok=True)

# -----------------------------------------
# FUNCIÓN GAUSSIANA
# -----------------------------------------
def gaussiana(x, A, mu, sigma, C):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + C


# -----------------------------------------
# FUNCIÓN PRINCIPAL
# -----------------------------------------
def medir_EW_gauss(filename):
    try:
        # -------------------------
        # LEER FITS (FULLDISKS)
        # -------------------------
        with fits.open(filename) as hdul:
            flux = hdul[2].data.field('FLUX')[0]
            wl_vac = hdul[2].data.field('WAVELENGTH')[0]

        wl = pyasl.vactoair2(wl_vac)

        # -------------------------
        # REGIÓN Hα (ajustada)
        # -------------------------
        mask = (wl >= 6555) & (wl <= 6570)
        wl_Ha = wl[mask]
        fl_Ha = flux[mask]

        if len(wl_Ha) < 10:
            return None

        # -------------------------
        # ESTIMACIONES INICIALES
        # -------------------------
        C_init = np.median(fl_Ha)
        A_init = np.max(fl_Ha) - C_init
        mu_init = wl_Ha[np.argmax(fl_Ha)]
        sigma_init = 1.0

        p0 = [A_init, mu_init, sigma_init, C_init]

        bounds = (
            [0, 6555, 0.2, 0],
            [np.inf, 6570, 10, np.inf]
        )

        popt, _ = curve_fit(
            gaussiana,
            wl_Ha,
            fl_Ha,
            p0=p0,
            bounds=bounds,
            maxfev=20000
        )

        A, mu, sigma, C = popt

        # -------------------------
        # EW EN λ (SIGNO CORRECTO)
        # -------------------------
        EW_lambda = - A * np.sqrt(2*np.pi) * sigma / C

        # -------------------------
        # VELOCIDAD
        # -------------------------
        vel = c * (wl_Ha - lambda_0) / lambda_0

        # -------------------------
        # EW EN VELOCIDAD
        # -------------------------
        EW_vel = (c / lambda_0) * EW_lambda

        nombre = os.path.basename(filename)

        print(f"{nombre}")
        print(f"   EW_lambda = {EW_lambda:.3f} Å")
        print(f"   EW_vel = {EW_vel:.3f} km/s")

        # =================================================
        # 📊 GRÁFICA COMBINADA (λ + velocidad)
        # =================================================
        plt.figure(figsize=(10,5))

        modelo = gaussiana(wl_Ha, *popt)

        # normalización visual
        plt.plot(wl_Ha, fl_Ha / C, 'k.', label="Datos")
        plt.plot(wl_Ha, modelo / C, 'r-', label="Ajuste")

        plt.axhline(1, linestyle='--', color='gray')

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

        plt.legend()
        plt.tight_layout()

        ruta = os.path.join(
            carpeta_graficas,
            nombre.replace(".fits", "_gauss.jpg")
        )

        plt.savefig(ruta, dpi=150)
        plt.close()

        return EW_lambda, EW_vel

    except Exception as e:
        print(f"{filename} -> ERROR: {e}")
        return None


# -----------------------------------------
# LOOP
# -----------------------------------------
results = []

carpeta = "FULLDISKS"

for file in os.listdir(carpeta):
    if file.endswith(".fits"):

        ruta = os.path.join(carpeta, file)

        res = medir_EW_gauss(ruta)

        if res is not None:
            EW_lambda, EW_vel = res

            results.append({
                "Archivo": file,
                "EW_lambda_A": EW_lambda,
                "EW_vel_kms": EW_vel
            })


# -----------------------------------------
# CSV
# -----------------------------------------
if results:
    df = pd.DataFrame(results)
    df.to_csv("EW_gauss.csv", index=False)
    print("\nCSV guardado: EW_gauss.csv")
else:
    print("\nNo hubo resultados válidos")

