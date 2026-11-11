import pandas as pd
import numpy as np

# ---------------------------------------------------------
# Constantes
# ---------------------------------------------------------

c = 3e5          # velocidad de la luz (km/s)
lambda_0 = 6562.8  # longitud de onda de Halpha (A)

# ---------------------------------------------------------
# Lectura de los resultados obtenidos previamente
# ---------------------------------------------------------

bandas = pd.read_csv("SoloY.csv")
gauss  = pd.read_csv("SoloY_EW_GAUSS.csv")

# ---------------------------------------------------------
# Seleccion de las tres primeras columnas
# ---------------------------------------------------------

bandas = bandas.iloc[:, :3].copy()
gauss  = gauss.iloc[:, :3].copy()

bandas.columns = ["Archivo", "EW_lambda_A", "Metodo"]
gauss.columns  = ["Archivo", "EW_lambda_A", "Metodo"]

# ---------------------------------------------------------
# Eliminacion de registros sin identificador o EW
# ---------------------------------------------------------

bandas = bandas.dropna(subset=["Archivo", "EW_lambda_A"])
gauss  = gauss.dropna(subset=["Archivo", "EW_lambda_A"])

# ---------------------------------------------------------
# Eliminacion de archivos repetidos
# Se conserva una unica medicion por espectro
# ---------------------------------------------------------

bandas = bandas.drop_duplicates(subset="Archivo")
gauss  = gauss.drop_duplicates(subset="Archivo")

# ---------------------------------------------------------
# Seleccion de los espectros presentes en ambos metodos
# ---------------------------------------------------------

comunes = set(bandas["Archivo"]) & set(gauss["Archivo"])

bandas = bandas[bandas["Archivo"].isin(comunes)].copy()
gauss  = gauss[gauss["Archivo"].isin(comunes)].copy()

# Ordenamiento para garantizar correspondencia
bandas = bandas.sort_values("Archivo").reset_index(drop=True)
gauss  = gauss.sort_values("Archivo").reset_index(drop=True)

# ---------------------------------------------------------
# Conversion de EW de A a km/s
#
# EW_v = (c / lambda_0) EW_lambda
# ---------------------------------------------------------

bandas["EW_vel_kms"] = (
    c / lambda_0
) * bandas["EW_lambda_A"]

gauss["EW_vel_kms"] = (
    c / lambda_0
) * gauss["EW_lambda_A"]

# ---------------------------------------------------------
# Renombrar columnas
# ---------------------------------------------------------

bandas = bandas.rename(columns={
    "EW_lambda_A": "EW_Bandas_A",
    "EW_vel_kms": "EW_Bandas_kms"
})

gauss = gauss.rename(columns={
    "EW_lambda_A": "EW_Gauss_A",
    "EW_vel_kms": "EW_Gauss_kms"
})

# ---------------------------------------------------------
# Union de los resultados
# ---------------------------------------------------------

resultado = pd.merge(
    bandas[["Archivo", "EW_Bandas_A", "EW_Bandas_kms"]],
    gauss[["Archivo", "EW_Gauss_A", "EW_Gauss_kms"]],
    on="Archivo",
    how="inner"
)

# ---------------------------------------------------------
# Diferencias entre los metodos
# ---------------------------------------------------------

resultado["Delta_EW_A"] = (
    resultado["EW_Bandas_A"]
    - resultado["EW_Gauss_A"]
)

resultado["Abs_Delta_EW_A"] = (
    resultado["Delta_EW_A"].abs()
)

resultado["Delta_EW_kms"] = (
    resultado["EW_Bandas_kms"]
    - resultado["EW_Gauss_kms"]
)

resultado["Abs_Delta_EW_kms"] = (
    resultado["Delta_EW_kms"].abs()
)

# ---------------------------------------------------------
# Guardar resultados
# ---------------------------------------------------------

resultado.to_csv(
    "EW_comparacion_lambda_velocidad.csv",
    index=False
)

print("Numero de espectros comunes:", len(resultado))
print("\nResultados guardados en:")
print("EW_comparacion_lambda_velocidad.csv")
