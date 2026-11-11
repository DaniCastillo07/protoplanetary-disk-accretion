# =====================================================
# MEDICIÓN AUTOMÁTICA DE EW(Hα)
# Daniela Castillo
# =====================================================

# -----------------------------------------------------
# 1️⃣ IMPORTAMOS LIBRERÍAS
# -----------------------------------------------------

import sys              # Para leer argumentos desde terminal
import os               # Para manejar carpetas
import numpy as np      # Matemáticas
import pandas as pd     # Crear tabla CSV
import matplotlib.pyplot as plt  # Graficar
from astropy.io import fits       # Abrir archivos FITS
from PyAstronomy import pyasl     # Convertir vacío → aire


# -----------------------------------------------------
# 2️⃣ VERIFICAR QUE EL USUARIO PASÓ EL .lis
# -----------------------------------------------------

if len(sys.argv) < 2:
    print("Uso correcto:")
    print("python3 ancho_equivalente.py fulldisk.lis")
    sys.exit()

lista_archivos = sys.argv[1]


# -----------------------------------------------------
# 3️⃣ CREAR CARPETAS DE SALIDA
# -----------------------------------------------------

carpeta_graficas = "GRAFICAS_EW"
carpeta_resultados = "RESULTADOS_EW"

os.makedirs(carpeta_graficas, exist_ok=True)
os.makedirs(carpeta_resultados, exist_ok=True)


# -----------------------------------------------------
# 4️⃣ LEER ARCHIVO .lis
# -----------------------------------------------------

if not os.path.isfile(lista_archivos):
    print("No existe:", lista_archivos)
    sys.exit()

with open(lista_archivos, "r") as f:
    archivos = f.read().splitlines()

# Asegurarnos que estén en FULLDISKS
archivos = [
    j if j.startswith("FULLDISKS/")
    else os.path.join("FULLDISKS", j)
    for j in archivos
]


# -----------------------------------------------------
# 5️⃣ PREPARAR TABLA DE RESULTADOS
# -----------------------------------------------------

nombres = []
ews = []


print("\n===== INICIANDO PROCESAMIENTO =====\n")


# -----------------------------------------------------
# 6️⃣ LOOP PRINCIPAL
# -----------------------------------------------------

for archivo in archivos:

    nombre_simple = os.path.basename(archivo)
    print("Procesando:", nombre_simple)

    try:

        # Abrimos FITS
        with fits.open(archivo) as hdul:
            fl = hdul[2].data.field('FLUX')[0]
            ivar = hdul[2].data.field('IVAR')[0]
            wl_vac = hdul[2].data.field('WAVELENGTH')[0]

        # Convertimos vacío → aire
        wl = pyasl.vactoair2(wl_vac, mode='ciddor')

        # -------------------------------------------------
        # DEFINIMOS BANDAS
        # -------------------------------------------------

        banda_linea = (6533, 6593)
        banda_cont_izq = (6533, 6545)
        banda_cont_der = (6580, 6593)

        # Seleccionamos regiones
        ind_linea = (wl >= banda_linea[0]) & (wl <= banda_linea[1])
        ind_izq = (wl >= banda_cont_izq[0]) & (wl <= banda_cont_izq[1])
        ind_der = (wl >= banda_cont_der[0]) & (wl <= banda_cont_der[1])

        wl_linea = wl[ind_linea]
        fl_linea = fl[ind_linea]

        if len(wl_linea) == 0:
            raise ValueError("No hay datos en Hα")

        # -------------------------------------------------
        # CONTINUO (PROMEDIO)
        # -------------------------------------------------

        cont_izq = np.mean(fl[ind_izq])
        cont_der = np.mean(fl[ind_der])
        continuo = (cont_izq + cont_der) / 2

        # Normalizamos
        fl_norm = fl_linea / continuo

        # -------------------------------------------------
        # CÁLCULO DE EW
        # -------------------------------------------------

        ew = np.trapz(1 - fl_norm, wl_linea)

        print(f"   EW = {ew:.3f} Å")

    except Exception as e:
        print("   Error:", e)
        ew = np.nan

    nombres.append(nombre_simple)
    ews.append(ew)


    # -------------------------------------------------
    # 7️⃣ GRAFICAR
    # -------------------------------------------------

    plt.figure(figsize=(10,5))

    # Espectro completo
    plt.plot(wl, fl/continuo, color='gray', alpha=0.4)

    # Línea Hα
    plt.plot(wl_linea, fl_norm, color='black')

    # Bandas sombreadas
    plt.axvspan(*banda_cont_izq, color='cyan', alpha=0.3)
    plt.axvspan(*banda_cont_der, color='cyan', alpha=0.3)
    plt.axvspan(*banda_linea, color='red', alpha=0.15)

    plt.axhline(1, linestyle='--')
    plt.axvline(6563, linestyle='--', color='blue')

    plt.title(f"{nombre_simple}   EW(Hα) = {ew:.2f} Å")
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flujo normalizado")

    plt.tight_layout()

    ruta_grafica = os.path.join(
        carpeta_graficas,
        nombre_simple.replace(".fits", "_EW.jpg")
    )

    plt.savefig(ruta_grafica, dpi=200)
    plt.close()


# -----------------------------------------------------
# 8️⃣ GUARDAR CSV
# -----------------------------------------------------

df = pd.DataFrame({
    "archivo": nombres,
    "EW_Halpha": ews
})

ruta_csv = os.path.join(
    carpeta_resultados,
    lista_archivos.replace(".lis", "_EW.csv")
)

df.to_csv(ruta_csv, index=False)

print("\n===== PROCESO TERMINADO =====")
print("Resultados en:", ruta_csv)
print("Gráficas en carpeta:", carpeta_graficas)

