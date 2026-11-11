import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ================================================================
# CONFIGURACION
# ================================================================

archivo_csv = "resultados_EW.csv"

carpeta_salida = "EW_resultados"

os.makedirs(carpeta_salida, exist_ok=True)


# ================================================================
# LECTURA DE LOS RESULTADOS
# ================================================================

df = pd.read_csv(archivo_csv)


# ================================================================
# IDENTIFICACION DE COLUMNAS
# ================================================================

col_band = "EW_lambda_A"
col_gauss = "EW_gauss_A"


# ================================================================
# SELECCION DE DATOS VALIDOS
# ================================================================

datos = df[[col_band, col_gauss]].dropna().copy()

EW_band = datos[col_band].to_numpy()
EW_gauss = datos[col_gauss].to_numpy()


# ================================================================
# CALCULO DE DIFERENCIAS
# ================================================================

datos["Delta_EW"] = (
    datos[col_band] - datos[col_gauss]
)

datos["Abs_Delta_EW"] = (
    np.abs(datos["Delta_EW"])
)


# ================================================================
# ESTADISTICOS DE LOS ANCHOS EQUIVALENTES
# ================================================================

estadisticos_EW = pd.DataFrame({
    "Estadistico": [
        "N",
        "Media",
        "Desviacion estandar",
        "Minimo",
        "Q1 (25 %)",
        "Mediana",
        "Q3 (75 %)",
        "Maximo"
    ],

    "$EW_{\\mathrm{Bandas}}$ (\\AA)": [
        len(EW_band),
        np.mean(EW_band),
        np.std(EW_band, ddof=1),
        np.min(EW_band),
        np.percentile(EW_band, 25),
        np.median(EW_band),
        np.percentile(EW_band, 75),
        np.max(EW_band)
    ],

    "$EW_{\\mathrm{Gauss}}$ (\\AA)": [
        len(EW_gauss),
        np.mean(EW_gauss),
        np.std(EW_gauss, ddof=1),
        np.min(EW_gauss),
        np.percentile(EW_gauss, 25),
        np.median(EW_gauss),
        np.percentile(EW_gauss, 75),
        np.max(EW_gauss)
    ]
})


# ================================================================
# ESTADISTICOS DE LAS DIFERENCIAS
# ================================================================

delta = datos["Delta_EW"].to_numpy()
abs_delta = datos["Abs_Delta_EW"].to_numpy()


estadisticos_delta = pd.DataFrame({
    "Estadistico": [
        "N",
        "Media",
        "Desviacion estandar",
        "Minimo",
        "Q1 (25 %)",
        "Mediana",
        "Q3 (75 %)",
        "Maximo"
    ],

    "$\\Delta EW$ (\\AA)": [
        len(delta),
        np.mean(delta),
        np.std(delta, ddof=1),
        np.min(delta),
        np.percentile(delta, 25),
        np.median(delta),
        np.percentile(delta, 75),
        np.max(delta)
    ],

    "$|\\Delta EW|$ (\\AA)": [
        len(abs_delta),
        np.mean(abs_delta),
        np.std(abs_delta, ddof=1),
        np.min(abs_delta),
        np.percentile(abs_delta, 25),
        np.median(abs_delta),
        np.percentile(abs_delta, 75),
        np.max(abs_delta)
    ]
})


# ================================================================
# DISTRIBUCION SEGUN LA DIFERENCIA ABSOLUTA
# ================================================================

umbrales = [0.5, 1, 2, 5, 10]

rangos = []

for limite in umbrales:

    cantidad = np.sum(abs_delta <= limite)

    porcentaje = (
        cantidad / len(abs_delta) * 100
    )

    rangos.append({
        "Criterio":
            rf"$|\Delta EW| \leq {limite}\,\AA$",

        "Numero de espectros":
            cantidad,

        "Porcentaje (%)":
            porcentaje
    })


tabla_rangos = pd.DataFrame(rangos)


# ================================================================
# TABLA COMPLETA DE LOS 208 ESPECTROS
# ================================================================

tabla_208 = datos.copy()

tabla_208.insert(
    0,
    "N",
    np.arange(1, len(tabla_208) + 1)
)

tabla_208.to_csv(
    os.path.join(
        carpeta_salida,
        "tabla_208_espectros.csv"
    ),
    index=False
)


# ================================================================
# GRAFICA 1:
# COMPARACION ENTRE METODOS
# ================================================================

plt.figure(figsize=(7, 7))

plt.scatter(
    EW_gauss,
    EW_band,
    s=25,
    alpha=0.7
)

limite_min = min(
    np.min(EW_gauss),
    np.min(EW_band)
)

limite_max = max(
    np.max(EW_gauss),
    np.max(EW_band)
)

plt.plot(
    [limite_min, limite_max],
    [limite_min, limite_max],
    "--",
    linewidth=1.5,
    label=r"$y=x$"
)

plt.xlabel(
    r"$EW_{\mathrm{Gauss}}$ (\AA)"
)

plt.ylabel(
    r"$EW_{\mathrm{Bandas}}$ (\AA)"
)

plt.title(
    r"Comparación de los métodos de medición de $EW$ de H$\alpha$"
)

plt.legend()

plt.grid(
    linestyle=":",
    alpha=0.5
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        carpeta_salida,
        "comparacion_EW_metodos.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ================================================================
# GRAFICA 2:
# DIFERENCIA ABSOLUTA
# ================================================================

plt.figure(figsize=(9, 5))

plt.scatter(
    np.arange(1, len(abs_delta) + 1),
    abs_delta,
    s=22,
    alpha=0.7
)

plt.axhline(
    np.median(abs_delta),
    linestyle="--",
    linewidth=1.5,
    label=(
        r"Mediana = "
        f"{np.median(abs_delta):.3f}"
        r" $\AA$"
    )
)

plt.xlabel(
    "Espectro"
)

plt.ylabel(
    r"$|\Delta EW|$ (\AA)"
)

plt.title(
    r"Diferencia absoluta entre los métodos de medición"
)

plt.legend()

plt.grid(
    linestyle=":",
    alpha=0.5
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        carpeta_salida,
        "diferencia_EW.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ================================================================
# EXPORTACION DE LAS TABLAS
# ================================================================

def formatear_latex(df_tabla):

    return df_tabla.to_latex(
        index=False,
        escape=False,
        float_format="%.3f"
    )


latex_EW = formatear_latex(
    estadisticos_EW
)

latex_delta = formatear_latex(
    estadisticos_delta
)

latex_rangos = formatear_latex(
    tabla_rangos
)


with open(
    os.path.join(
        carpeta_salida,
        "tablas_estadisticas_EW.tex"
    ),
    "w",
    encoding="utf-8"
) as archivo:

    archivo.write(
        "% Tabla de estadisticos de EW\n\n"
    )

    archivo.write(
        latex_EW
    )

    archivo.write(
        "\n\n% Tabla de diferencias\n\n"
    )

    archivo.write(
        latex_delta
    )

    archivo.write(
        "\n\n% Tabla de rangos\n\n"
    )

    archivo.write(
        latex_rangos
    )


# ================================================================
# MOSTRAR RESULTADOS
# ================================================================

print("\n===================================")
print("ESTADISTICOS DE EW")
print("===================================\n")

print(
    estadisticos_EW.to_string(
        index=False
    )
)


print("\n===================================")
print("ESTADISTICOS DE LAS DIFERENCIAS")
print("===================================\n")

print(
    estadisticos_delta.to_string(
        index=False
    )
)


print("\n===================================")
print("DISTRIBUCION DE LAS DIFERENCIAS")
print("===================================\n")

print(
    tabla_rangos.to_string(
        index=False
    )
)


print("\n===================================")
print(
    f"Analizados: {len(datos)} espectros"
)
print("===================================")
