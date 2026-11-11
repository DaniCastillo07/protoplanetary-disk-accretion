"""
Construye los 3 catalogos de entrada para MassAge (uno por formato de
temperatura disponible: SpT, lgTeff, Teff lineal), cruzando cada fuente de
Teff con la fotometria del catalogo base.

Requiere: teff_consolidado.csv (o los 4 archivos fuente), y
KINCANDS_2MASS_W4_or_M24_fulldisk.cat

Genera:
- massage_cat1_SpT.csv       (46 estrellas, optspec)
- massage_cat2_lgTeff.csv    (46 estrellas, apogee)
- massage_cat3_Teff.csv      (57 estrellas, GBM + TIC v8.2)
"""
import pandas as pd
import numpy as np

DATA_DIR = './'

def loadsrc(fn):
    with open(DATA_DIR + fn) as f:
        header_line = f.readline()
    cols = header_line.lstrip('#').split()
    return pd.read_csv(DATA_DIR + fn, sep=r'\s+', skiprows=1, header=None, names=cols)

opt = loadsrc('fulldisk_optspec.txt')
gbm = loadsrc('fulldisk_GBM.txt')
apo = loadsrc('fulldisk_apogee.txt')
tic = loadsrc('fulldisk_TICv82.txt')

with open(DATA_DIR + 'KINCANDS_2MASS_W4_or_M24_fulldisk.cat') as f:
    header_line = f.readline()
bcols = header_line.lstrip('#').split()
base = pd.read_csv(DATA_DIR + 'KINCANDS_2MASS_W4_or_M24_fulldisk.cat', sep=r'\s+',
                    skiprows=1, header=None, names=bcols, quotechar='"')


def photo_block(gaia_ids):
    """Extrae y renombra las columnas de fotometria/astrometria requeridas
    por MassAge para un conjunto de gaia_source_id dado."""
    sub = base[base['source_id_1'].isin(gaia_ids)].copy()
    return pd.DataFrame({
        'Source': sub['source_id_1'],
        'RA': sub['ra_epoch2000'],
        'DEC': sub['dec_epoch2000'],
        'BPmag': sub['Bpmag'],
        'Gmag': sub['Gmag'],
        'RPmag': sub['Rpmag'],
        'Jmag': sub['Jmag'],
        'Hmag': sub['Hmag'],
        'e_BPmag': sub['Bperr'],
        'e_Gmag': sub['Gerr'],
        'e_RPmag': sub['Rperr'],
        'e_Jmag': sub['e_Jmag'],
        'e_Hmag': sub['e_Hmag'],
        'PLXcorr': sub['PLXcorr'],
        'e_Plx': sub['PLXerr'],
    })


# --- Catalogo 1: optspec -> SpT / SpTerr ---
p1 = photo_block(opt['gaia_source_id'])
t1 = opt[['gaia_source_id', 'SpT', 'SpTe']].rename(
    columns={'gaia_source_id': 'Source', 'SpTe': 'SpTerr'})
cat1 = p1.merge(t1, on='Source')
assert len(cat1) == 46
cat1.to_csv('massage_cat1_SpT.csv', index=False)

# --- Catalogo 2: apogee -> lgTeff / elgTeff ---
p2 = photo_block(apo['gaia_source_id'])
t2 = apo[['gaia_source_id', 'logTeff', 'e_logTeff']].rename(
    columns={'gaia_source_id': 'Source', 'logTeff': 'lgTeff', 'e_logTeff': 'elgTeff'})
cat2 = p2.merge(t2, on='Source')
assert len(cat2) == 46
cat2.to_csv('massage_cat2_lgTeff.csv', index=False)

# --- Catalogo 3: GBM + TIC v8.2 -> Teff / eTeff (lineal) ---
p3a = photo_block(gbm['gaia_source_id'])
t3a = gbm[['gaia_source_id', 'Teff_GBM', 'eTeff_GBM']].rename(
    columns={'gaia_source_id': 'Source', 'Teff_GBM': 'Teff', 'eTeff_GBM': 'eTeff'})
cat3a = p3a.merge(t3a, on='Source')

p3b = photo_block(tic['gaia_source_id'])
t3b = tic[['gaia_source_id', 'Teff', 's_Teff']].rename(
    columns={'gaia_source_id': 'Source', 's_Teff': 'eTeff'})
cat3b = p3b.merge(t3b, on='Source')

cat3 = pd.concat([cat3a, cat3b], ignore_index=True)
assert len(cat3) == 57
cat3.to_csv('massage_cat3_Teff.csv', index=False)

print('cat1 (SpT):', cat1.shape)
print('cat2 (lgTeff):', cat2.shape)
print('cat3 (Teff):', cat3.shape)
print('Total:', len(cat1) + len(cat2) + len(cat3))
