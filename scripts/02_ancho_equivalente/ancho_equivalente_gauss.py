# =====================================================
# MEDICIÓN DE EW(Hα) POR AJUSTE GAUSSIANO (EMISIÓN)
# =====================================================

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.optimize import curve_fit


# -----------------------------------------------------
# FUNCIÓN GAUSSIANA (continuo normalizado = 1)
# -----------------------------------------------------

def gaussiana(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + 1


# -----------------------------------------------------
# VERIFICAR ARGUMENTO
# -----------------------------------------------------

if len(sys.argv) < 2:
    print("Uso:")
    print("python3 ancho_equivalente_gauss.py lista.csv")
    sys.exit()

archivo_lista = sys.argv[1]


# -----------------------------------------------------
# CARPETAS
# -----------------------------------------------------

carpeta_graficas = "GRAFICAS_EW_GAUSS"
carpeta_resultados = "RESULTADOS_EW_GAUSS"

os.makedirs(carpeta_graficas, exist_ok=True)
os.makedirs(carpeta_resultados, exist_ok=True)


# -----------------------------------------------------
# LEER LISTA
# -----------------------------------------------------

tabla = pd.read_csv(archivo_lista, header=None)
lista_espectros = tabla[0]


# -----------------------------------------------------
# RESULTADOS
# -----------------------------------------------------

nombres = []
ews = []

print("\n===== AJUSTE GAUSSIANO Hα =====\n")


# =====================================================
# LOOP PRINCIPAL
# =====================================================

for nombre in lista_espectros:

    ruta = nombre.strip()

    print("\nProcesando:", ruta)

    try:

        # -------------------------------------------------
        # LEER FITS
        # -------------------------------------------------

        with fits.open(ruta) as hdul:
            fl = hdul[2].data.field("FLUX")[0]
            wl_vac = hdul[2].data.field("WAVELENGTH")[0]

        wl = pyasl.vactoair2(wl_vac, mode="ciddor")


        # -------------------------------------------------
        # CONTINUO (MISMO DEL MÉTODO DE BANDAS)
        # -------------------------------------------------

        mask_left = (wl >= 6490) & (wl <= 6510)
        mask_right = (wl >= 6600) & (wl <= 6620)

        if np.sum(mask_left) == 0 or np.sum(mask_right) == 0:
            raise ValueError("No hay puntos suficientes para estimar el continuo.")

        cont_left = np.median(fl[mask_left])
        cont_right = np.median(fl[mask_right])

        continuo = np.median([cont_left, cont_right])

        if continuo <= 0 or not np.isfinite(continuo):
            raise ValueError("Continuo inválido.")

        fl_norm = fl / continuo


        # -------------------------------------------------
        # REGIÓN Hα
        # -------------------------------------------------

        mask = (wl >= 6530) & (wl <= 6595)

        wl_Ha = wl[mask]
        fl_Ha = fl_norm[mask]

        if len(wl_Ha) < 10:
            raise ValueError("Muy pocos puntos.")


        print(f"Continuo = {continuo:.4f}")
        print(f"Flujo máximo normalizado = {np.max(fl_Ha):.2f}")


        # -------------------------------------------------
        # PARÁMETROS INICIALES
        # -------------------------------------------------

        A_init = np.max(fl_Ha) - 1

        mu_init = wl_Ha[np.argmax(fl_Ha)]

        sigma_init = 1.5

        p0 = [A_init, mu_init, sigma_init]


        # -------------------------------------------------
        # AJUSTE
        # -------------------------------------------------

        bounds = (
            [0, 6550, 0.1],
            [100, 6575, 10]
        )

        popt, pcov = curve_fit(
            gaussiana,
            wl_Ha,
            fl_Ha,
            p0=p0,
            bounds=bounds,
            maxfev=30000
        )

        A, mu, sigma = popt

        print(f"A = {A:.3f}")
        print(f"mu = {mu:.3f}")
        print(f"sigma = {sigma:.3f}")


        # -------------------------------------------------
        # EW
        # -------------------------------------------------

        EW = -A * np.sqrt(2*np.pi) * sigma


        # -------------------------------------------------
        # CONTROL DE CALIDAD
        # -------------------------------------------------

        if sigma > 8:
            raise ValueError(f"Sigma anómala: {sigma:.2f}")

        if abs(EW) > 300:
            raise ValueError(f"EW anómalo: {EW:.2f}")


        modelo = gaussiana(wl_Ha, *popt)

        print(f"EW = {EW:.3f} Å")


        nombres.append(os.path.basename(ruta))
        ews.append(EW)


        # -------------------------------------------------
        # GRÁFICA
        # -------------------------------------------------

        plt.figure(figsize=(10,5))

        plt.plot(
            wl_Ha,
            fl_Ha,
            "k.",
            ms=4,
            label="Datos normalizados"
        )

        plt.plot(
            wl_Ha,
            modelo,
            "r-",
            lw=2,
            label="Ajuste gaussiano"
        )

        plt.axhline(
            1,
            color="gray",
            ls="--"
        )

        plt.axvline(
            6562.8,
            color="blue",
            ls="--",
            label="Hα"
        )

        plt.title(
            f"{os.path.basename(ruta)}\nEW = {EW:.2f} Å"
        )

        plt.xlabel("Longitud de onda (Å)")
        plt.ylabel("Flujo normalizado")

        plt.legend()

        plt.tight_layout()

        nombre_fig = os.path.basename(ruta).replace(
            ".fits",
            "_GAUSS.jpg"
        )

        plt.savefig(
            os.path.join(carpeta_graficas, nombre_fig),
            dpi=200
        )

        plt.close()


    except Exception as e:

        print("Error:", e)
        continue


# -----------------------------------------------------
# GUARDAR CSV
# -----------------------------------------------------

df = pd.DataFrame({

    "archivo": nombres,

    "EW_Halpha_Gauss": ews

})

ruta_csv = os.path.join(

    carpeta_resultados,

    archivo_lista.replace(".csv", "_EW_GAUSS.csv")

)

df.to_csv(ruta_csv, index=False)

print("\n===== TERMINADO =====")

print("CSV:", ruta_csv)

print("Gráficas:", carpeta_graficas)
