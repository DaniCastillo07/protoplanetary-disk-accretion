# Códigos de la tesis — Acreción de discos protoplanetarios en Orión

Repositorio de scripts de Python utilizados para el procesamiento y análisis de espectros LAMOST y la caracterización de estrellas T Tauri en la región de Orión. Corresponde al trabajo de tesis de **Daniela Castillo Valderrama**.

## Organización

| Carpeta | Contenido | Apéndices relacionados |
| --- | --- | --- |
| `scripts/01_inspeccion_espectral/` | Inspección visual y clasificación manual de espectros | A |
| `scripts/02_ancho_equivalente/` | EW de Hα por bandas y ajuste gaussiano; conversión al espacio de velocidades; comparación de métodos e histogramas | B, C, E, G, H |
| `scripts/03_ancho_W10/` | Determinación de la anchura de Hα al 10 % de su máximo | D, I |
| `scripts/04_diagrama_HR/` | Diagrama H-R de evolución pre-secuencia principal basado en valores tabulados usados en la tesis | F |
| `scripts/05_parametros_estelares/` | Consolidación de temperatura efectiva; preparación de catálogos de entrada y consolidación de resultados de MassAge | J y análisis complementarios |
| `scripts/06_tasa_acrecion/` | Estimación de la tasa de acreción y comparación con parámetros estelares y exceso infrarrojo | K, L |

## Scripts recuperados de los apéndices

Estos cuatro archivos `.py` se extrajeron del texto LaTeX de los apéndices proporcionado para esta colección:

- `E_conversion_EW_velocidades.py` — conversión conjunta de EW y cálculo de diferencias entre métodos.
- `F_diagrama_HR_Baraffe.py` — visualización de trayectorias e isócronas.
- `G_comparacion_metodos_EW.py` — estadísticas, gráficos y tablas de comparación de EW.
- `H_histograma_EW.py` — histogramas superpuestos de EW por los dos métodos.

**Modificación mínima:** en `G_comparacion_metodos_EW.py` se unió una cadena de texto de un título gráfico que estaba dividida incorrectamente entre dos líneas, para resolver un error de sintaxis. Los cálculos no se modificaron.

Los demás scripts originales se conservaron con sus nombres y contenido.

## Dependencias

Para instalar las bibliotecas de Python utilizadas por los scripts:

```bash
python -m pip install -r requirements.txt
```

**Datos no incluidos:** los scripts necesitan archivos FITS y tablas o catálogos de entrada (`.csv`, `.cat`, `.txt`, `.lis`, `.head`) que no forman parte de este repositorio. Algunas rutas, nombres de archivos y nombres de columnas son específicos del entorno original; deben revisarse antes de ejecutar cada etapa. Las entradas no están unificadas automáticamente entre todos los scripts.

**MassAge:** el repositorio incluye únicamente scripts del análisis de la tesis que generan entradas o procesan salidas de MassAge; **no redistribuye el código interno de la herramienta MassAge**.

