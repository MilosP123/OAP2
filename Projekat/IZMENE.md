# Primenjene izmene notebook-a — 20. septembar 2026.

Polazni GitHub commit: `64f1807ff55681d23a1e0188af11a2bb5023dc93`. Brojevi su fizičke pozicije ćelija (1–33), uključujući Markdown. nbqa svoje kodne ćelije može numerisati drugačije. Sačuvani stari izlazi su uklonjeni.

Ranija lokalna revizija od 17. septembra je pregledana i dopunjena. Ovo poređenje prikazuje original sa GitHub-a i sadašnji kod ove revizije.

## Ćelija 1

Dodat uvod: pitanje, opseg, nepromenljivost originalnih podataka i status rezultata.

### Original
```python

```

### Izmenjeno (markdown)
```markdown
# Okruženje galaksija i AGN aktivnost — TNG100-1

Pitanje: kako se AGN aktivnost razlikuje sa masom matičnog halo-a, uz kontrolu
zvezdane mase i razdvajanje centralnih i satelitskih galaksija?
Analiza koristi snapshot 99 (z ≈ 0), bez analize prečki.

Pokrenuti **Restart & Run All** iz direktorijuma `Projekat`. Potreban je ceo
stvarni group katalog u `data/group99/`; instalacija paketa je opisana u README.
Notebook ne zamenjuje nedostajuće podatke simuliranim. Sintetički podaci koriste
se isključivo u odvojenim testovima i njihovi izlazi se ne objavljuju kao naučni.

Stari sačuvani izlazi uklonjeni su. Novi naučni rezultati postoje tek nakon
uspešnog izvršavanja nad stvarnim podacima. HDF5 ulazi otvaraju se samo za čitanje.

```

## Ćelija 2

Svi importi su u prvoj kodnoj ćeliji; uklonjeni neupotrebljeni importi, dodati alati za poreklo izlaza.

### Original
```python
import numpy as np
import h5py as hp # biblioteka za manipulaciiju nad hdf fajlovima
import pandas as pd
import glob
import os
import illustris_python as il
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rc
from scipy.stats import pearsonr, spearmanr
```

### Izmenjeno (code)
```python
"""Reproduktivna analiza okruženja i AGN aktivnosti u TNG100-1."""
import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import h5py as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr, spearmanr

```

## Ćelija 3

Relativne putanje, snapshot 99 i fiksan seed bootstrap-a.

### Original
```python
GROUPCAT_DIR = "/Volumes/Tina/Milos/TNG100/group99"
```

### Izmenjeno (code)
```python
GROUPCAT_DIR = Path("data/group99")
OUTPUT_DIR = Path("outputs")
SNAP = 99
RANDOM_SEED = 42
N_BOOTSTRAP = 1000

```

## Ćelija 4

Numerički delovi, kompletnost, konzistentna zaglavlja i h/a/z/BoxSize; jasne poruke za nedostajuće podatke.

### Original
```python
files = sorted(glob.glob(os.path.join(GROUPCAT_DIR, "*.hdf5")))
with hp.File(files[0], "r") as data:
    print("Glavne grupe u fajlu:")
    print(list(data.keys()))

    print("\nHeader atributi:")
    for key in data["Header"].attrs.keys():
        print(key, "=", data["Header"].attrs[key])

    print("\nPolja u Group:")
    print(list(data["Group"].keys())[:30])

    print("\nPolja u Subhalo:")
    print(list(data["Subhalo"].keys())[:30])
```

### Izmenjeno (code)
```python
# Čitamo samo delove izabranog snapshota, ne proizvoljne HDF5 fajlove.
files = list(GROUPCAT_DIR.glob(f"groupcat-{SNAP}.*.hdf5"))
files += list(GROUPCAT_DIR.glob(f"fof_subhalo_tab_{SNAP:03d}.*.hdf5"))
if not files:
    raise FileNotFoundError(
        f"Nema delova kataloga u {GROUPCAT_DIR}. Potrebni su svi delovi; vidi README.md."
    )
try:
    files.sort(key=lambda name: int(name.name.split(".")[-2]))
except ValueError as exc:
    raise ValueError("Broj dela HDF5 kataloga mora biti ceo broj u nazivu fajla.") from exc

count_keys = ("NumFiles", "Ngroups_Total", "Nsubgroups_Total",
              "Ngroups_ThisFile", "Nsubgroups_ThisFile")
physical_keys = ("HubbleParam", "Redshift", "Time", "BoxSize")
header = None
input_parts = []
for filename in files:
    with hp.File(filename, "r") as data:
        if "Header" not in data:
            raise KeyError(f"Nedostaje Header: {filename}")
        attrs = dict(data["Header"].attrs)
        for key in count_keys + physical_keys:
            if key not in attrs:
                raise KeyError(f"Nedostaje Header/{key}: {filename}")
            value = np.asarray(attrs[key])
            if value.ndim != 0 or value.dtype.kind not in "iuf":
                raise ValueError(f"Header/{key} mora biti numerički skalar: {filename}")
            if not np.isfinite(value):
                raise ValueError(f"Nevažeći Header/{key}: {filename}")
        for key in count_keys:
            if attrs[key] < 0 or attrs[key] != int(attrs[key]):
                raise ValueError(f"Header/{key} mora biti nenegativan ceo broj: {filename}")
        if header is None:
            header = attrs
        for key in count_keys[:3]:
            if attrs[key] != header[key]:
                raise ValueError(f"Neusklađen Header/{key}: {filename}")
        for key in physical_keys:
            if not np.isclose(attrs[key], header[key], rtol=1e-10, atol=1e-12):
                raise ValueError(f"Pomešani katalozi, Header/{key}: {filename}")
        input_parts.append({
            "name": filename.name, "size_bytes": filename.stat().st_size,
            "groups": int(attrs["Ngroups_ThisFile"]),
            "subhalos": int(attrs["Nsubgroups_ThisFile"]),
        })

expected_files = int(header["NumFiles"])
expected_groups = int(header["Ngroups_Total"])
expected_subhalos = int(header["Nsubgroups_Total"])
chunk_ids = [int(name.name.split(".")[-2]) for name in files]
if chunk_ids != list(range(expected_files)):
    raise ValueError("Nedostaju delovi kataloga ili postoje duplikati brojeva delova.")
if expected_groups == 0 or expected_subhalos == 0:
    raise ValueError("Katalog nema halo grupe ili subhaloe za analizu.")
if sum(part["groups"] for part in input_parts) != expected_groups:
    raise ValueError("Zbir Ngroups_ThisFile ne odgovara Ngroups_Total.")
if sum(part["subhalos"] for part in input_parts) != expected_subhalos:
    raise ValueError("Zbir Nsubgroups_ThisFile ne odgovara Nsubgroups_Total.")
if SNAP != 99 or not np.isclose(header["Redshift"], 0., atol=1e-8):
    raise ValueError("Ovaj notebook je namenjen snapshot-u 99, z približno 0.")
if not np.isclose(header["Time"], 1 / (1 + header["Redshift"]), atol=1e-8):
    raise ValueError("Header/Time nije saglasan sa Redshift.")
if (not np.isclose(header["HubbleParam"], 0.6774)
        or not np.isclose(header["BoxSize"], 75000.)):
    raise ValueError("h ili BoxSize ne odgovaraju očekivanom TNG100 katalogu.")
print("Broj delova:", len(files))
print("Očekivani broj grupa/subhaloa:", expected_groups, expected_subhalos)

```

## Ćelija 5

Pragovi zadržani uz eksplicitno operativno značenje.

### Original
```python
SNAP = 99
# Pragovi za podelu okruženja
FIELD_GROUP_LIMIT = 5e12       # Msun
GROUP_CLUSTER_LIMIT = 8e13     # Msun, motivisano Paul et al. 2017
# Prag za minimalnu zvezdanu masu galaksije
MSTAR_MIN = 1e9                # Msun
# Prag za AGN definiciju
LAMBDA_EDD_LIMIT = 1e-3
```

### Izmenjeno (code)
```python
FIELD_GROUP_LIMIT = 5e12       # Msun; operativna donja granica grupa
GROUP_CLUSTER_LIMIT = 8e13     # Msun; motivisano Paul et al. (2017)
MSTAR_MIN = 1e9                # Msun; proveriti osetljivost na ovaj izbor
LAMBDA_EDD_LIMIT = 1e-3        # operativna definicija AGN kandidata

```

## Ćelija 6

Objašnjene ulazne tabele i jedinice u Markdown-u.

### Original
```python
#Treba da napravim histogram podataka
#Za tvoj projekat Subhalo je mesto gde se nalaze galaksije, dok je Group mesto gde se nalaze njihova okruženja
#Statisticki raspodela edingtonovog odnosa i mase crne rupe za z=0,1,2 (podaci za akrecioni disk0) Eta i luminoznost proveriti
```

### Izmenjeno (markdown)
```markdown
## Učitavanje i jedinice

`Group` sadrži halo grupe, `Subhalo` njihove subhaloe. Brojevi redova su globalni
indeksi samo kada su svi delovi spojeni numeričkim redosledom.
Proveravamo dimenzije, ukupan broj objekata i konzistentnost veza.

Mase su u `1e10 Msun/h`; akrecija je u `(1e10 Msun/h)/(0.978 Gyr/h)`.
`SubhaloMassType[:, 4]` daje ukupnu vezanu zvezdanu masu. BH masa i akrecija
su sume po subhalou. SFR je već u `Msun/yr`, položaji u `ckpc/h`.
Izvor: [TNG specifikacije](https://www.tng-project.org/data/docs/specifications/).

```

## Ćelija 7

Provera punih dimenzija i celobrojnih indeksa; nedostajuće neprazno polje je greška; unsigned sentinel -1; bezbedan logaritam.

### Original
```python
#Ovo je zapravo kako cemo da podatke iz svakog fajla spajamo u jedan niz, pored toga proverava da li postoje podaci
def ucitaj_polje(files, grupa, polje):

    lista = []

    for file in files:
        with hp.File(file, "r") as data:
            if grupa in data and polje in data[grupa]:
                lista.append(data[grupa][polje][:])

    if len(lista) == 0:
        raise ValueError(f"Nisam našao polje {grupa}/{polje}")

    return np.concatenate(lista, axis=0)
```

### Izmenjeno (code)
```python
def ucitaj_polje(fajlovi, grupa, polje):
    """Spoji polje uz provere dimenzija, tipa i broja redova u svakom delu."""
    keys = {"Group": ("Ngroups_ThisFile", "Ngroups_Total"),
            "Subhalo": ("Nsubgroups_ThisFile", "Nsubgroups_Total")}
    count_key, total_key = keys[grupa]
    tail_shape = {"SubhaloMassType": (6,), "SubhaloPos": (3,)}.get(polje, ())
    integer_fields = {"SubhaloGrNr", "GroupFirstSub", "GroupNsubs", "SubhaloFlag"}
    chunks = []
    expected_total = None
    for part in fajlovi:
        with hp.File(part, "r") as handle:
            count = int(handle["Header"].attrs[count_key])
            total = int(handle["Header"].attrs[total_key])
            if expected_total is None:
                expected_total = total
            elif total != expected_total:
                raise ValueError(f"Neusklađen ukupan broj: {part}")
            if grupa not in handle or polje not in handle[grupa]:
                if count == 0:
                    continue
                raise KeyError(f"Nedostaje {grupa}/{polje} u {part}")
            dataset = handle[grupa][polje]
            expected_shape = (count,) + tail_shape
            if not isinstance(dataset, hp.Dataset) or dataset.shape != expected_shape:
                raise ValueError(f"Pogrešan oblik {grupa}/{polje} u {part}; "
                                 f"očekivano {expected_shape}")
            if dataset.dtype.kind not in "iuf":
                raise ValueError(f"Polje {grupa}/{polje} nije numeričko: {part}")
            if polje in integer_fields and dataset.dtype.kind not in "iu":
                raise ValueError(f"Polje {grupa}/{polje} mora biti celobrojno: {part}")
            values = dataset[:]
            # Dokumentacija dozvoljava unsigned zapis sentinela -1 za GroupFirstSub.
            if polje == "GroupFirstSub" and values.dtype.kind == "u":
                sentinel = values == np.iinfo(values.dtype).max
                if np.any(values[~sentinel] > np.iinfo(np.int64).max):
                    raise ValueError("GroupFirstSub indeks izvan int64 opsega.")
                values = values.astype(np.int64)
                values[sentinel] = -1
            elif polje in integer_fields:
                if values.dtype.kind == "u" and np.any(values > np.iinfo(np.int64).max):
                    raise ValueError(f"{polje} indeks izvan int64 opsega.")
                values = values.astype(np.int64)
            chunks.append(values)
    if not chunks:
        raise ValueError(f"Nisam našao polje {grupa}/{polje}")
    result = np.concatenate(chunks, axis=0)
    if len(result) != expected_total:
        raise ValueError(f"Nepotpun katalog za {grupa}/{polje}")
    return result


def pozitivan_log(values):
    """Log10 samo za pozitivne konačne vrednosti, ostalo je NaN."""
    values = np.asarray(values, dtype=np.float64)
    valid = np.isfinite(values) & (values > 0)
    return np.log10(values, where=valid, out=np.full(values.shape, np.nan))

```

## Ćelija 8

Korišćenje već validiranog zaglavlja.

### Original
```python
#Ovde sada za neki faj ucitavamo pojedinacne vrednosti a svi ovi podaci inace treba da budu isti
with hp.File(files[0], "r") as data:
    h = data["Header"].attrs["HubbleParam"]
    redshift = data["Header"].attrs["Redshift"]
    scale_factor = data["Header"].attrs["Time"]

print("h =", h)
print("redshift =", redshift)
print("scale factor =", scale_factor)
```

### Izmenjeno (code)
```python
h = float(header["HubbleParam"])
redshift = float(header["Redshift"])
scale_factor = float(header["Time"])
print("h, z, a =", h, redshift, scale_factor)

```

## Ćelija 9

Zadržana tri Group polja; sažeti tačni komentari.

### Original
```python
#Ucitavanje podataka za grupu ("Group")
Group_M_Crit200 = ucitaj_polje(files, "Group", "Group_M_Crit200") #Preko nje cemo kasnije da odradimo klasifikaciju (field,group,cluster)
#Group_M_Crit200 je unutar radijusa gde je srednja gustina 200 puta veća od kritične gustine Univerzuma
GroupFirstSub = ucitaj_polje(files, "Group", "GroupFirstSub") #Ucitava indeks prvog subhaloa u svakoj grupi (obicno centralna galaksija), moze da nam posluzi za odredjivanje satelit galaksija
GroupNsubs = ucitaj_polje(files, "Group", "GroupNsubs") #Ucitava broj subhalo objekata

```

### Izmenjeno (code)
```python
# Masa u sferi srednje gustine 200 * kritična gustina.
Group_M_Crit200 = ucitaj_polje(files, "Group", "Group_M_Crit200")
GroupFirstSub = ucitaj_polje(files, "Group", "GroupFirstSub")
GroupNsubs = ucitaj_polje(files, "Group", "GroupNsubs")

```

## Ćelija 10

Zadržana Subhalo polja, dodat domen SubhaloFlag 0/1.

### Original
```python
#Ucitavanje podataka za subhalo ("Subhalo")
SubhaloFlag = ucitaj_polje(files, "Subhalo", "SubhaloFlag") #Ovo ucitava flag koji govori da li je subhalo pogodan za standardnu analizu
SubhaloGrNr = ucitaj_polje(files, "Subhalo", "SubhaloGrNr") #SubhaloGrNr za svaki subhalo kaže u kojoj Group se nalazi
SubhaloMassType = ucitaj_polje(files, "Subhalo", "SubhaloMassType") #Ovo učitava masu subhalo-a podeljenu po tipovima čestica. Dvodimenzioni niz kod koga nam je bitna kolona 4 koja prica o masi cestica
SubhaloBHMass = ucitaj_polje(files, "Subhalo", "SubhaloBHMass") #Ovo učitava ukupnu masu crnih rupa u svakom subhalou (potrebno za edingtonov odnos)
SubhaloBHMdot = ucitaj_polje(files, "Subhalo", "SubhaloBHMdot") #Ovo je akreciona stopa crne rupe, tj. koliko materije crna rupa akretuje po jedinici vremena (isto za edingtonov odnos)
SubhaloSFR = ucitaj_polje(files, "Subhalo", "SubhaloSFR") #Ovo učitava star formation rate, odnosno stopu formiranja zvezda u galaksiji (nije toliko vazno trenutno za nas)
SubhaloPos = ucitaj_polje(files, "Subhalo", "SubhaloPos") #Ovo učitava pozicije subhalo-a u simulacionoj kutiji (nije neophodno)
```

### Izmenjeno (code)
```python
SubhaloFlag = ucitaj_polje(files, "Subhalo", "SubhaloFlag")
SubhaloGrNr = ucitaj_polje(files, "Subhalo", "SubhaloGrNr")
SubhaloMassType = ucitaj_polje(files, "Subhalo", "SubhaloMassType")
SubhaloBHMass = ucitaj_polje(files, "Subhalo", "SubhaloBHMass")
SubhaloBHMdot = ucitaj_polje(files, "Subhalo", "SubhaloBHMdot")
SubhaloSFR = ucitaj_polje(files, "Subhalo", "SubhaloSFR")
SubhaloPos = ucitaj_polje(files, "Subhalo", "SubhaloPos")
if not np.isin(SubhaloFlag, [0, 1]).all():
    raise ValueError("SubhaloFlag mora sadržati samo 0 ili 1.")

```

## Ćelija 11

Jedinstven naziv M200c_all i float64 pre konverzije; akrecija bez dodatnog h.

### Original
```python
# Masa halo-a u Msun
M200c_a = Group_M_Crit200 * 1e10 / h

# Zvezdana masa galaksije
Mstar = SubhaloMassType[:, 4] * 1e10 / h

# Masa crne rupe
MBH = SubhaloBHMass * 1e10 / h

#Akrecionu stopa crne rupe 
Mdot_BH = SubhaloBHMdot * 1e10 / 0.978e9
```

### Izmenjeno (code)
```python
# Float64 pre množenja. h se poništava samo u konverziji akrecione stope.
M200c_all = Group_M_Crit200.astype(np.float64) * 1e10 / h
Mstar = SubhaloMassType[:, 4].astype(np.float64) * 1e10 / h
MBH = SubhaloBHMass.astype(np.float64) * 1e10 / h
Mdot_BH = SubhaloBHMdot.astype(np.float64) * 1e10 / 0.978e9

```

## Ćelija 12

Jedinstven M200c_parent, opseg indeksa, broj članova, globalni redosled i početak svake grupe; centrale/sateliti pre rezova.

### Original
```python
#Za svaku galaksiju pronalazimo masu halo-a u kom se ona nalazi
parent_group_id = SubhaloGrNr
SubfindID = np.arange(len(SubhaloGrNr))

valid_parent = parent_group_id < len(M200c_a)

M200c_p = np.full(len(SubhaloGrNr), np.nan)
M200c_p[valid_parent] = M200c_all[parent_group_id[valid_parent]]
```

### Izmenjeno (code)
```python
parent_group_id = SubhaloGrNr
SubfindID = np.arange(len(SubhaloGrNr), dtype=np.int64)
valid_parent = (parent_group_id >= 0) & (parent_group_id < len(M200c_all))
if not valid_parent.all():
    raise ValueError("Nevažeći SubhaloGrNr u kompletnom katalogu.")
if np.any(GroupNsubs < 0):
    raise ValueError("GroupNsubs mora biti nenegativan.")
counts = np.bincount(parent_group_id, minlength=len(GroupNsubs))
if not np.array_equal(counts, GroupNsubs):
    raise ValueError("SubhaloGrNr i GroupNsubs se ne slažu.")
if np.any(np.diff(parent_group_id) < 0):
    raise ValueError("SubhaloGrNr nije u globalnom redosledu halo grupa.")

occupied_groups = np.flatnonzero(GroupNsubs > 0)
central_ids = GroupFirstSub[occupied_groups]
if np.any((central_ids < 0) | (central_ids >= len(SubfindID))):
    raise ValueError("Nevažeći GroupFirstSub za nepraznu grupu.")
if np.any(GroupFirstSub[GroupNsubs == 0] != -1):
    raise ValueError("Prazna grupa mora imati GroupFirstSub = -1.")
expected_first = np.cumsum(GroupNsubs, dtype=np.int64) - GroupNsubs
if not np.array_equal(central_ids, expected_first[occupied_groups]):
    raise ValueError("GroupFirstSub nije prvi globalni subhalo svoje grupe.")
if not np.array_equal(parent_group_id[central_ids], occupied_groups):
    raise ValueError("GroupFirstSub i SubhaloGrNr se ne slažu.")

M200c_parent = M200c_all[parent_group_id]
is_central = SubfindID == GroupFirstSub[parent_group_id]
is_satellite = ~is_central
if int(is_central.sum()) != len(occupied_groups):
    raise ValueError("Svaka neprazna grupa mora imati tačno jednu centralnu galaksiju.")

```

## Ćelija 13

Konačne fizički validne BH vrednosti; nulta akrecija ostaje nula; overflow ostaje nepoznat; ETA=0.1 zadržano i dokumentovano.

### Original
```python
MSUN_G = 1.98847e33
YEAR_S = 365.25 * 24 * 3600
C_LIGHT = 2.99792458e10

#Pretvaramo u decimalni zapis
MBH = MBH.astype(np.float64)
Mdot_BH = Mdot_BH.astype(np.float64)

#Radijativna efikasnost akrecije
ETA = 0.1

# Mdot iz Msun/yr u g/s
Mdot_g_s = Mdot_BH * MSUN_G / YEAR_S

# Bolometrijska luminoznost
L_bol = ETA * Mdot_g_s * C_LIGHT**2

# Eddingtonova luminoznost
L_edd = 1.26e38 * MBH

# Eddingtonov odnos
lambda_edd = np.full_like(MBH, np.nan, dtype=float)

#Uslov za koji se dodeljuju vrednosti (Ako objekat ima crnu rupu, ali trenutno ne akretuje, onda je: Mdot_BH = 0)
mask_bh = MBH > 0
lambda_edd[mask_bh] = L_bol[mask_bh] / L_edd[mask_bh]
```

### Izmenjeno (code)
```python
MSUN_G = 1.98847e33
YEAR_S = 365.25 * 24 * 3600
C_LIGHT = 2.99792458e10
ETA = 0.1  # Zadržana postprocesna pretpostavka; TNG model koristi 0.2.
if not (0 < ETA <= 1 and LAMBDA_EDD_LIMIT > 0 and MSTAR_MIN > 0
        and 0 < FIELD_GROUP_LIMIT < GROUP_CLUSTER_LIMIT):
    raise ValueError("Nevažeći pragovi ili radijativna efikasnost.")

# Ako nevažeći ekstrem izazove overflow, vrednost ostaje nepoznata (ne AGN).
with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
    Mdot_g_s = Mdot_BH * MSUN_G / YEAR_S
    L_bol = ETA * Mdot_g_s * C_LIGHT**2
    L_edd = 1.26e38 * MBH
    lambda_edd = np.full_like(MBH, np.nan, dtype=float)
    mask_bh = (np.isfinite(MBH) & (MBH > 0)
               & np.isfinite(Mdot_BH) & (Mdot_BH >= 0)
               & np.isfinite(L_bol) & np.isfinite(L_edd))
    lambda_edd[mask_bh] = L_bol[mask_bh] / L_edd[mask_bh]
lambda_edd[~np.isfinite(lambda_edd)] = np.nan
# MBH > 0 i Mdot_BH = 0 daju lambda_edd = 0: validna neaktivnost.

```

## Ćelija 14

Formule, uzorak i ograničenja naučnih definicija.

### Original
```python

```

### Izmenjeno (markdown)
```markdown
## Definicije uzorka i aktivnosti

$L_{bol}=\eta\dot M c^2$, $L_{Edd}=1.26\times10^{38}(M_{BH}/M_\odot)$ erg/s,
$\lambda=L_{bol}/L_{Edd}$. `ETA=0.1` je eksplicitna postprocesna pretpostavka,
ne merena efikasnost niti automatski stvarna luminoznost pri svakoj akreciji.
TNG model koristi 0.2; tabelom osetljivosti proveravamo oba izbora.

Osnovni uzorak: SubhaloFlag=1, konačna Mstar ≥ 1e9 Msun, konačna M200c > 0.
Rezovi po BH aktivnosti ne ulaze u izbor osnovnog uzorka. Nulta akrecija ostaje
u imeniocu frakcije; nepoznat status se broji odvojeno.
`field` je ovde operativna masena kategorija, ne potvrda prostorne izolovanosti.
Granica group/cluster-like motivisana je radom
[Paul et al. (2017)](https://arxiv.org/pdf/1706.01916).

```

## Ćelija 15

Nevažeća halo masa dobija unknown; tačne granice kategorija.

### Original
```python
#Sada pravimo klasifikaciju okruzenja
okruzenje = np.full(len(M200c_p), "", dtype=object)

okruzenje[M200c_p < FIELD_GROUP_LIMIT] = "field"

okruzenje[(M200c_p >= FIELD_GROUP_LIMIT) & (M200c_p < GROUP_CLUSTER_LIMIT)] = "group"

okruzenje[M200c_p >= GROUP_CLUSTER_LIMIT] = "cluster-like"
```

### Izmenjeno (code)
```python
okruzenje = np.full(len(M200c_parent), "unknown", dtype=object)
valid_halo_mass = np.isfinite(M200c_parent) & (M200c_parent > 0)
okruzenje[valid_halo_mass & (M200c_parent < FIELD_GROUP_LIMIT)] = "field"
okruzenje[valid_halo_mass & (M200c_parent >= FIELD_GROUP_LIMIT)
         & (M200c_parent < GROUP_CLUSTER_LIMIT)] = "group"
okruzenje[valid_halo_mass & (M200c_parent >= GROUP_CLUSTER_LIMIT)] = "cluster-like"

```

## Ćelija 16

Nullable AGN status, čisti logaritmi, clean bez selekcije po aktivnosti; tabela kumulativnih rezova.

### Original
```python
#Pravimo tabelu

catalog = pd.DataFrame()
catalog["SubfindID"] = SubfindID
catalog["parent_group_id"] = parent_group_id
catalog["SubhaloFlag"] = SubhaloFlag
catalog["Mstar"] = Mstar
catalog["logMstar"] = np.log10(Mstar, where=(Mstar > 0), out=np.full_like(Mstar, np.nan))
catalog["M200c_parent"] = M200c_p
catalog["logM200c_parent"] = np.log10(M200c_p,where=(M200c_p > 0),out=np.full_like(M200c_p, np.nan))
catalog["MBH"] = MBH
catalog["Mdot_BH"] = Mdot_BH
catalog["lambda_edd"] = lambda_edd
catalog["SFR"] = SubhaloSFR
catalog["environment"] = okruzenje
catalog["x"] = SubhaloPos[:, 0]
catalog["y"] = SubhaloPos[:, 1]
catalog["z"] = SubhaloPos[:, 2]
#logaritamske vrednosti
catalog["logMBH"] = np.nan
catalog.loc[catalog["MBH"] > 0, "logMBH"] = np.log10(catalog.loc[catalog["MBH"] > 0, "MBH"])
catalog["logMdot_BH"] = np.nan
catalog.loc[catalog["Mdot_BH"] > 0, "logMdot_BH"] = np.log10(catalog.loc[catalog["Mdot_BH"] > 0, "Mdot_BH"])
catalog["log_lambda_edd"] = np.nan
catalog.loc[catalog["lambda_edd"] > 0, "log_lambda_edd"] = np.log10(catalog.loc[catalog["lambda_edd"] > 0, "lambda_edd"])
#Da li je agn kandidat
catalog["is_agn"] = catalog["lambda_edd"] >= LAMBDA_EDD_LIMIT
#LAMBDA_EDD_LIMIT = 1e-4
#LAMBDA_EDD_LIMIT = 1e-3
#LAMBDA_EDD_LIMIT = 1e-2
#Cisnenje kataloga
clean = catalog[(catalog["SubhaloFlag"] == 1) &(catalog["Mstar"] >= MSTAR_MIN) &(catalog["M200c_parent"] > 0) &np.isfinite(catalog["M200c_parent"])].copy()

```

### Izmenjeno (code)
```python
catalog = pd.DataFrame({
    "SubfindID": SubfindID, "parent_group_id": parent_group_id,
    "SubhaloFlag": SubhaloFlag, "is_central": is_central,
    "is_satellite": is_satellite, "Mstar": Mstar, "M200c_parent": M200c_parent,
    "MBH": MBH, "Mdot_BH": Mdot_BH, "lambda_edd": lambda_edd,
    "SFR": SubhaloSFR, "environment": okruzenje,
    "x": SubhaloPos[:, 0], "y": SubhaloPos[:, 1], "z": SubhaloPos[:, 2],
})
for column, log_column in (("Mstar", "logMstar"), ("M200c_parent", "logM200c_parent"),
                           ("MBH", "logMBH"), ("Mdot_BH", "logMdot_BH"),
                           ("lambda_edd", "log_lambda_edd")):
    catalog[log_column] = pozitivan_log(catalog[column])

catalog["is_agn"] = pd.Series(pd.NA, index=catalog.index, dtype="boolean")
known_bh = np.isfinite(lambda_edd)
catalog.loc[known_bh, "is_agn"] = lambda_edd[known_bh] >= LAMBDA_EDD_LIMIT
no_bh = np.isfinite(MBH) & (MBH == 0) & np.isfinite(Mdot_BH) & (Mdot_BH == 0)
catalog.loc[no_bh, "is_agn"] = False

# Kumulativni izveštaj o selekciji, uključujući nepoznato okruženje.
base_mask = np.ones(len(catalog), dtype=bool)
selection_rows = []
selection_steps = [
    ("svi_subhaloi", np.ones(len(catalog), dtype=bool)),
    ("SubhaloFlag_1", SubhaloFlag == 1),
    ("konacna_Mstar", np.isfinite(Mstar)),
    ("Mstar_min", Mstar >= MSTAR_MIN),
    ("validna_M200c", valid_halo_mass),
]
for step, condition in selection_steps:
    base_mask &= condition
    row = {"korak": step, "ukupno": int(base_mask.sum())}
    for env in ("field", "group", "cluster-like", "unknown"):
        row[env] = int(np.sum(base_mask & (okruzenje == env)))
    selection_rows.append(row)
selection = pd.DataFrame(selection_rows)
clean = catalog.loc[base_mask].copy()
print(selection.to_string(index=False))
print("Nepoznat AGN status u osnovnom uzorku:", clean["is_agn"].isna().sum())
# x, y, z su ckpc/h; fizički kpc = koordinata * a / h.

```

## Ćelija 17

Poznati/ukupni imenilac, broj halo-a, granice zbog nepoznatih; centrale/sateliti i uzorak sa BH.

### Original
```python
#Ovaj kod pravi tabelu koja za svako okruženje pokazuje koliko ima galaksija, koliko ima AGN-a, kolika je AGN frakcija,
#i koje su tipične vrednosti zvezdane mase, mase halo-a, Eddingtonovog odnosa i BH akrecione stope
summary = clean.groupby("environment").agg(
    broj_galaksija=("SubfindID", "count"),
    broj_AGN=("is_agn", "sum"),
    AGN_frakcija=("is_agn", "mean"),
    medijana_logMstar=("logMstar", "median"),
    medijana_logM200c=("logM200c_parent", "median"),
    medijana_log_lambda_edd=("log_lambda_edd", "median"),
    medijana_logMdot_BH=("logMdot_BH", "median")
)
print(summary)

#lambda_edd ne predstavlja medijanu svih galaksija, nego praktično medijanu onih koje imaju definisanu pozitivnu vrednost
```

### Izmenjeno (code)
```python
def pregled_uzorka(frame, by):
    """Frakcije sa eksplicitnim imeniocem i granicama za nepoznate statuse."""
    result = frame.groupby(by, observed=True).agg(
        broj_galaksija=("SubfindID", "count"),
        broj_haloa=("parent_group_id", "nunique"),
        broj_sa_AGN_statusom=("is_agn", "count"),
        broj_AGN=("is_agn", "sum"),
        AGN_frakcija=("is_agn", "mean"),
        medijana_logMstar=("logMstar", "median"),
        medijana_logM200c=("logM200c_parent", "median"),
        medijana_log_lambda_edd=("log_lambda_edd", "median"),
        medijana_logMdot_BH=("logMdot_BH", "median"),
    )
    result["broj_nepoznatih"] = result["broj_galaksija"] - result["broj_sa_AGN_statusom"]
    result["AGN_frakcija_min"] = result["broj_AGN"] / result["broj_galaksija"]
    result["AGN_frakcija_max"] = (
        (result["broj_AGN"] + result["broj_nepoznatih"]) / result["broj_galaksija"]
    )
    return result


summary = pregled_uzorka(clean, "environment")
summary_central = pregled_uzorka(clean, ["environment", "is_central"])
summary_bh = pregled_uzorka(clean[np.isfinite(clean["MBH"]) & (clean["MBH"] > 0)],
                          "environment")
print(summary)
print(summary_central)
print("Frakcije su među poznatim statusima; min/max opisuju nepoznate, ne statistički CI.")
print("Medijane log(lambda) i log(Mdot) uključuju samo pozitivne konačne vrednosti.")

```

## Ćelija 18

Konačni parovi, premali/konstantni nizovi, vraćanje korelacija za CSV i objašnjenje p=0.

### Original
```python
#Primena Pirona i Spermana kao statosticka procena korelacija
def korelacija(df, x, y):
    pom = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()

    print("\n===================================")
    print("Korelacija:", x, "vs", y)
    print("Broj tačaka:", len(pom))

    if len(pom) < 3:
        print("Nema dovoljno podataka.")
        return

    pearson_r, pearson_p = pearsonr(pom[x], pom[y])
    spearman_r, spearman_p = spearmanr(pom[x], pom[y])

    print("Pearson r =", pearson_r)
    print("Pearson p =", pearson_p)

    print("Spearman rho =", spearman_r)
    print("Spearman p =", spearman_p)
```

### Izmenjeno (code)
```python
def korelacija(frame, x_col, y_col):
    """Deskriptivne korelacije konačnih parova; naivni p nije halo-aware test."""
    pairs = frame[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    result = {"x": x_col, "y": y_col, "n": len(pairs),
              "pearson_r": np.nan, "pearson_p": np.nan,
              "spearman_rho": np.nan, "spearman_p": np.nan}
    if len(pairs) < 3 or pairs[x_col].nunique() < 2 or pairs[y_col].nunique() < 2:
        print("Korelacija nije definisana: premalo parova ili konstantan niz.")
        return result
    result["pearson_r"], result["pearson_p"] = pearsonr(pairs[x_col], pairs[y_col])
    result["spearman_rho"], result["spearman_p"] = spearmanr(pairs[x_col], pairs[y_col])
    print(result)
    if result["pearson_p"] == 0 or result["spearman_p"] == 0:
        print("p=0 je numerički ispis, ne dokaz apsolutne sigurnosti.")
    return result

```

## Ćelija 19

Korelacije se čuvaju kao tabela.

### Original
```python
korelacija(clean, "logM200c_parent", "log_lambda_edd")
korelacija(clean, "logM200c_parent", "logMdot_BH")
korelacija(clean, "logMstar", "log_lambda_edd")
```

### Izmenjeno (code)
```python
correlations = pd.DataFrame([
    korelacija(clean, "logM200c_parent", "log_lambda_edd"),
    korelacija(clean, "logM200c_parent", "logMdot_BH"),
    korelacija(clean, "logMstar", "log_lambda_edd"),
])

```

## Ćelija 20

Ograničenja korelacija i interpretacija grafikona.

### Original
```python

```

### Izmenjeno (markdown)
```markdown
## Grafikoni i statistička interpretacija

Korelacije opisuju povezanost, ne uzročnost. Njihove uobičajene p-vrednosti ne
uvažavaju zajednički halo, razlike u Mstar i udelu satelita.
Logaritamski grafikoni prikazuju samo pozitivne konačne parove, uključujući
objekte ispod AGN praga. Nulte akrecije su i dalje uključene u AGN frakcije.

Za poređenje sirovog i odabranog uzorka koristimo iste ivice binova.
Prazan skup pozitivnih parova daje objašnjenje u grafikonu, ne grešku.

```

## Ćelija 21

AGN prag, stvarni broj parova, prazni grafik i zatvaranje figure; relativni izlaz.

### Original
```python
plot_data = clean.dropna(subset=["logM200c_parent", "log_lambda_edd"])

plt.figure(figsize=(7, 5))

plt.scatter(plot_data["logM200c_parent"],plot_data["log_lambda_edd"],s=8,alpha=0.5)
plt.axvline(np.log10(FIELD_GROUP_LIMIT), linestyle="--", label="field/group")
plt.axvline(np.log10(GROUP_CLUSTER_LIMIT), linestyle="-.", label="group/cluster-like")
plt.xlabel(r"$\log_{10}(M_{200c}/M_\odot)$")
plt.ylabel(r"$\log_{10}(\lambda_{\rm Edd})$")
plt.title("AGN aktivnost u zavisnosti od mase okruženja")

plt.legend()
plt.tight_layout()
plt.savefig("lambda_edd_vs_m200c.png", dpi=200)
plt.show()
```

### Izmenjeno (code)
```python
plot_data = clean[["logM200c_parent", "log_lambda_edd"]].dropna()
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(plot_data["logM200c_parent"], plot_data["log_lambda_edd"], s=8, alpha=0.5)
ax.axvline(np.log10(FIELD_GROUP_LIMIT), linestyle="--", label="field/group")
ax.axvline(np.log10(GROUP_CLUSTER_LIMIT), linestyle="-.", label="group/cluster-like")
ax.axhline(np.log10(LAMBDA_EDD_LIMIT), color="black", linestyle=":", label="AGN prag")
ax.set_xlabel(r"$\log_{10}(M_{200c}/M_\odot)$")
ax.set_ylabel(r"$\log_{10}(\lambda_{\rm Edd})$")
ax.set_title(f"Pozitivna akrecija u osnovnom uzorku: N={len(plot_data)}")
if plot_data.empty:
    ax.text(0.5, 0.5, "Nema pozitivnih konačnih parova", transform=ax.transAxes, ha="center")
ax.legend()
fig.tight_layout()
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT_DIR / "lambda_edd_vs_m200c.png", dpi=200)
plt.show()
plt.close(fig)

```

## Ćelija 22

Zajednička funkcija za postojeće 2D histograme, zaštita praznog LogNorm grafikona.

### Original
```python
plot_data = clean.dropna(subset=["logMBH", "log_lambda_edd"])

plt.figure(figsize=(7, 5))

plt.hist2d(plot_data["logMBH"],plot_data["log_lambda_edd"],bins=60)

plt.xlabel(r"$\log_{10}(M_{\rm BH}/M_\odot)$")
plt.ylabel(r"$\log_{10}(\lambda_{\rm Edd})$")
plt.title("2D raspodela: Eddingtonov odnos vs masa crne rupe")

plt.colorbar(label="Broj galaksija")

plt.tight_layout()
plt.show()
```

### Izmenjeno (code)
```python
def histogram_bh(frame, bins, title, output_name, logarithmic=False):
    """2D histogram sa istim granicama i bez LogNorm greške za prazan uzorak."""
    pairs = frame[["logMBH", "log_lambda_edd"]].dropna()
    figure, axis = plt.subplots(figsize=(7, 5))
    if pairs.empty:
        axis.text(0.5, 0.5, "Nema pozitivnih konačnih parova",
                  transform=axis.transAxes, ha="center")
    else:
        norm = LogNorm(vmin=1) if logarithmic else None
        histogram = axis.hist2d(pairs["logMBH"], pairs["log_lambda_edd"],
                                bins=bins, norm=norm, cmin=1 if logarithmic else None)
        figure.colorbar(histogram[3], ax=axis, label="Broj subhaloa")
    axis.set_xlabel(r"$\log_{10}(M_{\rm BH}/M_\odot)$")
    axis.set_ylabel(r"$\log_{10}(\lambda_{\rm Edd})$")
    axis.set_title(f"{title}, N={len(pairs)}")
    figure.tight_layout()
    figure.savefig(OUTPUT_DIR / output_name, dpi=200)
    plt.show()
    plt.close(figure)

```

## Ćelija 23

Iste ivice binova za poređenje sirovih i čistih podataka.

### Original
```python
#import numpy as np
#import matplotlib.pyplot as plt

#def plot_hist(data, name, bins=50, logx=False, logy=False):
#    data = np.asarray(data).ravel()
#    data = data[np.isfinite(data)]

#    if logx:
#        data = data[data > 0]
#        data = np.log10(data)
#        xlabel = f"log10({name})"
#    else:
#        xlabel = name

#    plt.figure(figsize=(7, 4))
#    plt.hist(data, bins=bins, edgecolor="black", alpha=0.7)

#    if logy:
#        plt.yscale("log")

#   plt.xlabel(xlabel)
#    plt.ylabel("Broj objekata")
#    plt.title(f"Histogram: {name}")
#    plt.tight_layout()
#    plt.show()
```

### Izmenjeno (code)
```python
raw_pairs = catalog[["logMBH", "log_lambda_edd"]].dropna()
if raw_pairs.empty:
    bh_bins = [np.linspace(0., 1., 51), np.linspace(-8., 0., 51)]
else:
    bh_bins = [np.histogram_bin_edges(raw_pairs[column], bins=50)
               for column in ("logMBH", "log_lambda_edd")]

```

## Ćelija 24

Sirovi 2D histogram sa jasnim naslovom i posebnim fajlom.

### Original
```python
#gas_mass = SubhaloMassType[:, 0]
#dm_mass = SubhaloMassType[:, 1]
#stellar_mass = SubhaloMassType[:, 4]
#bh_mass = SubhaloMassType[:, 5]


#plot_hist(gas_mass, "Gas Mass", logx=True, logy=True)
#plot_hist(dm_mass, "Dark Matter Mass", logx=True, logy=True)
#plot_hist(stellar_mass, "Stellar Mass", logx=True, logy=True)
#plot_hist(bh_mass, "Black Hole Mass", logx=True, logy=True)

#plot_hist(Mdot_g_s,"Mdot" , logx=True, logy=True)
#plot_hist(lambda_edd[mask_bh],"Edingtonov odnos",logx=True, logy=True)
#plot_hist(L_bol,"Bolimetrijska luminoznost",logx=True, logy=True)
#plot_hist(L_edd,"Eddingtonova luminoznost",logx=True, logy=True)
```

### Izmenjeno (code)
```python
histogram_bh(catalog, bh_bins, "Svi subhaloi, pozitivna akrecija", "bh_raw.png")

```

## Ćelija 25

Čisti 2D histogram sa istim binovima.

### Original
```python
#import numpy as np
#import matplotlib.pyplot as plt


x = np.array(MBH)
y = np.array(lambda_edd)

mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

plt.hist2d(np.log10(x[mask]), np.log10(y[mask]), bins=50)

plt.xlabel("log10(MBH)")
plt.ylabel("log10(Edingtonov odnos)")
plt.colorbar(label="Broj tačaka")
plt.show()
```

### Izmenjeno (code)
```python
histogram_bh(clean, bh_bins, "Osnovni uzorak, pozitivna akrecija", "bh_clean.png")

```

## Ćelija 26

LogNorm prikaz istog čistog uzorka; bez duplog grafikona.

### Original
```python
#ocisceni podaci
plot_data = clean.dropna(subset=["logMBH", "log_lambda_edd"])

plt.figure(figsize=(7, 5))

plt.hist2d(plot_data["logMBH"],plot_data["log_lambda_edd"],bins=50)

plt.xlabel(r"$\log_{10}(M_{\rm BH}/M_\odot)$")
plt.ylabel(r"$\log_{10}(\lambda_{\rm Edd})$")
plt.colorbar(label="Broj galaksija")

plt.tight_layout()
plt.show()
```

### Izmenjeno (code)
```python
histogram_bh(clean, bh_bins, "Osnovni uzorak, logaritamska skala brojanja",
             "bh_clean_lognorm.png", logarithmic=True)

```

## Ćelija 27

Frakcije prema Mstar i tipu galaksije; bootstrap celih halo-a; obuhvaćena i maksimalna masa na ivici bina.

### Original
```python
plot_data = clean.dropna(subset=["logMBH", "log_lambda_edd"])

plt.figure(figsize=(7, 5))

plt.hist2d(plot_data["logMBH"],plot_data["log_lambda_edd"],bins=50)

plt.xlabel(r"$\log_{10}(M_{\rm BH}/M_\odot)$")
plt.ylabel(r"$\log_{10}(\lambda_{\rm Edd})$")
plt.colorbar(label="Broj galaksija")

plt.tight_layout()
plt.show()
```

### Izmenjeno (code)
```python
def halo_interval(frame, rng, n_resamples=1000):
    """Percentilni 95% CI za odnos brojeva, resampliranjem celih halo grupa."""
    known = frame.dropna(subset=["is_agn"])
    groups = known.groupby("parent_group_id")["is_agn"].agg(["sum", "count"])
    if len(groups) < 2:
        return np.nan, np.nan, len(groups)
    counts_array = groups[["sum", "count"]].to_numpy(dtype=float)
    fractions = np.empty(n_resamples)
    for index in range(n_resamples):
        draw = rng.integers(0, len(groups), size=len(groups))
        numerator, denominator = counts_array[draw].sum(axis=0)
        fractions[index] = numerator / denominator
    lo, hi = np.quantile(fractions, [0.025, 0.975])
    return lo, hi, len(groups)


rng = np.random.default_rng(RANDOM_SEED)
mass_rows = []
if not clean.empty:
    mass_edges = np.arange(np.floor(clean["logMstar"].min() * 2) / 2,
                           np.ceil(clean["logMstar"].max() * 2) / 2 + 1., .5)
    mass_sample = clean.assign(mass_bin=pd.cut(clean["logMstar"], mass_edges, right=False))
    for (env, central, mass_bin), subset in mass_sample.groupby(
            ["environment", "is_central", "mass_bin"], observed=True):
        known_count = int(subset["is_agn"].count())
        agn_count = int(subset["is_agn"].sum())
        lo, hi, halo_count = halo_interval(subset, rng, N_BOOTSTRAP)
        mass_rows.append({
            "environment": env, "is_central": bool(central),
            "logMstar_mid": mass_bin.mid, "N": len(subset),
            "N_known": known_count, "N_AGN": agn_count, "N_halo_known": halo_count,
            "AGN_fraction": agn_count / known_count if known_count else np.nan,
            "ci_low": lo, "ci_high": hi,
        })
mass_fractions = pd.DataFrame(mass_rows, columns=[
    "environment", "is_central", "logMstar_mid", "N", "N_known", "N_AGN",
    "N_halo_known", "AGN_fraction", "ci_low", "ci_high",
])
print(mass_fractions)
print("Bootstrap CI je eksplorativan pri malom broju halo-a; ne meri kosmičku varijansu.")

```

## Ćelija 28

Ključni graf po masi i okruženju, odvojeno centrale/sateliti, sa 95% halo bootstrap intervalima.

### Original
```python
#Ovo nije dobro
#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm

#x = np.array(Mdot_g_s)
#y = np.array(lambda_edd)

#mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

#plt.hist2d(
#    np.log10(x[mask]),
#    np.log10(y[mask]),
#    bins=50,
#    norm=LogNorm()
#)

#plt.xlabel("log10(Mdot)")
#plt.ylabel("log10(Edingtonov odnos)")
#plt.colorbar(label="Broj tačaka")
#plt.show()
```

### Izmenjeno (code)
```python
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, central, label in zip(axes, [True, False], ["Centralne", "Satelitske"]):
    for env in ("field", "group", "cluster-like"):
        rows = mass_fractions[(mass_fractions["environment"] == env)
                              & (mass_fractions["is_central"] == central)]
        if not rows.empty:
            line, = ax.plot(rows["logMstar_mid"], rows["AGN_fraction"], "o-", label=env)
            # CI ne mora biti simetričan niti centriran na tačkastu procenu.
            ax.vlines(rows["logMstar_mid"], rows["ci_low"], rows["ci_high"],
                      color=line.get_color())
    ax.set_title(label)
    ax.set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
    ax.set_ylim(-0.03, 1.03)
    if ax.lines:
        ax.legend()
    else:
        ax.text(.5, .5, "Nema galaksija u uzorku", transform=ax.transAxes, ha="center")
axes[0].set_ylabel("AGN frakcija među poznatim statusima")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "agn_fraction_vs_mstar.png", dpi=200)
plt.show()
plt.close(fig)

```

## Ćelija 29

Tabela osetljivosti na ETA, AGN prag i minimalnu zvezdanu masu.

### Original
```python
#ocisceni podaci
x = np.array(clean["MBH"])
y = np.array(clean["lambda_edd"])

mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

plt.figure(figsize=(7, 5))

plt.hist2d(
    np.log10(x[mask]),
    np.log10(y[mask]),
    bins=50,
    norm=LogNorm()
)

plt.xlabel("log10(MBH)")
plt.ylabel("log10(Eddingtonov odnos)")
plt.colorbar(label="Broj tačaka")
plt.show()
```

### Izmenjeno (code)
```python
# Osetljivost na eksplicitne postprocesne pretpostavke, bez promene clean uzorka.
sensitivity_rows = []
for eta_test in (0.1, 0.2):
    for mass_cut in (1e9, 10**9.5, 1e10):
        for threshold in (1e-4, 1e-3, 1e-2):
            for env, subset in clean[clean["Mstar"] >= mass_cut].groupby("environment"):
                known = subset["is_agn"].notna()
                scaled_lambda = subset["lambda_edd"] * eta_test / ETA
                n_known = int(known.sum())
                n_agn = int((known & (scaled_lambda >= threshold)).sum())
                sensitivity_rows.append({
                    "eta": eta_test, "Mstar_min": mass_cut, "lambda_limit": threshold,
                    "environment": env, "N": len(subset), "N_known": n_known,
                    "N_AGN": n_agn, "AGN_fraction": n_agn / n_known if n_known else np.nan,
                })
sensitivity = pd.DataFrame(sensitivity_rows, columns=[
    "eta", "Mstar_min", "lambda_limit", "environment", "N", "N_known", "N_AGN",
    "AGN_fraction",
])

```

## Ćelija 30

Šta intervali i poređenja podržavaju, a koja ograničenja ostaju.

### Original
```python
#Neocisceni podaci
#x = np.array(MBH)
#y = np.array(lambda_edd)

#mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

#plt.hist2d(
#    np.log10(x[mask]),
#    np.log10(y[mask]),
#    bins=50,
#    norm=LogNorm()
#)

#plt.xlabel("log10(MBH)")
#plt.ylabel("log10(Edingtonov odnos)")
#plt.colorbar(label="Broj tačaka")
#plt.show()

```

### Izmenjeno (markdown)
```markdown
## Kako tumačiti rezultate

`summary` daje osnovnu demografiju; `summary_central` razdvaja tipove galaksija,
`summary_bh` uslovljava uzorak prisustvom rezolvirane BH. Imenilac svake frakcije
je eksplicitno naveden. Granice min/max za nepoznate statuse nisu intervali
poverenja. `mass_fractions` poredi mase po binovima od 0.5 dex; to smanjuje,
ali ne uklanja potpuno, razlike masenih raspodela unutar binova.

Bootstrap bira cele halo-e sa ponavljanjem. Grupe sa mnogo satelita zadržavaju
svoje težine u frakciji galaksija. CI sa <2 halo-a nije definisan; čak i uz
nekoliko halo-a može biti nestabilan. Prostorne korelacije između halo-a,
rezolucija, BH seeding i konačna zapremina ostaju ograničenja.

Tabela `sensitivity` omogućava proveru da li znak i veličina razlike zavise od
AGN praga, ETA ili zvezdane mase. Ne bira se prag naknadno zato što daje
najmanju p-vrednost. Nema zaključka o uzročnosti samo na osnovu korelacije.

```

## Ćelija 31

Reprodukcija i eksplicitno razdvajanje naučnih od sintetičkih test podataka.

### Original
```python
#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm

#x = np.array(MBH)
#y = np.array(Mdot_g_s)

#mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

#plt.hist2d(
#    np.log10(x[mask]),
#    np.log10(y[mask]),
#    bins=60,
#    norm=LogNorm()
#)

#plt.xlabel("log10(MBH)")
#plt.ylabel("log10(Mdot)")
#plt.colorbar(label="Broj tačaka")
#plt.show()
```

### Izmenjeno (markdown)
```markdown
## Reprodukcija i predaja

Izlazi sledeće ćelije pripadaju isključivo ovom izvršavanju. CSV osnovnog uzorka
zadržava i neaktivne galaksije; nije dovoljan za smanjenje masenog reza ispod
MSTAR_MIN. Izvoz sadrži kontrolnu sumu i opis metoda.

Za ispit: potrebni su stvarni podaci uz repo ili automatsko preuzimanje, uspešan
`pytest --nbmake --nbmake-timeout=60 Untitled7.ipynb` i `nbqa pylint Untitled7.ipynb`.
Testovi nad sintetičkim podacima proveravaju softver, ne predstavljaju zamenu za
stvarni uzorak. README opisuje šta još zahteva podatke sa eksternog diska.

```

## Ćelija 32

Kontrolni kriterijumi za smislenost rezultata.

### Original
```python
#x = np.array(Mdot_g_s)
#y = np.array(MBH)

#mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

#plt.hist2d(
#    np.log10(x[mask]),
#    np.log10(y[mask]),
#    bins=60,
#    norm=LogNorm()
#)

#plt.xlabel("log10(Mdot)")
#plt.ylabel("log10(MBH)")
#plt.colorbar(label="Broj tačaka")
#plt.show()
```

### Izmenjeno (markdown)
```markdown
## Provera smislenosti

- Zbir brojeva po okruženjima mora odgovarati osnovnom uzorku.
- Broj AGN ≤ broj poznatih statusa ≤ broj galaksija; frakcije su između 0 i 1.
- Nulta akrecija nije nedostajući podatak. Prazan graf pozitivnih vrednosti
  može biti smislen i mora biti jasno označen.
- U punom katalogu svaka neprazna grupa ima jednu centralnu galaksiju;
  u filtriranom uzorku centrala može biti izbačena bez promene identiteta satelita.
- Stare brojeve iz istorijskog pregleda ne koristiti kao očekivane rezultate testa.

```

## Ćelija 33

Izvoz svih tabela, podataka i metapodataka: jedinice, pragovi, delovi ulaza, verzije, seed, hash izvedenog kataloga.

### Original
```python

```

### Izmenjeno (code)
```python
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
catalog_path = OUTPUT_DIR / "catalog_base.csv.gz"
clean.to_csv(catalog_path, index=False, compression={"method": "gzip", "mtime": 0})
for name, table in (("summary", summary), ("summary_central", summary_central),
                    ("summary_bh", summary_bh)):
    table.to_csv(OUTPUT_DIR / f"{name}.csv")
for name, table in (("selection", selection), ("correlations", correlations),
                    ("mass_fractions", mass_fractions), ("sensitivity", sensitivity)):
    table.to_csv(OUTPUT_DIR / f"{name}.csv", index=False)
metadata = {
    "schema_version": 1, "analysis_revision": "2026-09-20",
    "simulation_expected": "TNG100-1", "snapshot": SNAP,
    "h": h, "redshift": redshift, "scale_factor": scale_factor,
    "eta": ETA, "mstar_min_msun": MSTAR_MIN, "lambda_edd_limit": LAMBDA_EDD_LIMIT,
    "field_group_limit_msun": FIELD_GROUP_LIMIT,
    "group_cluster_limit_msun": GROUP_CLUSTER_LIMIT,
    "seed": RANDOM_SEED, "bootstrap_resamples": N_BOOTSTRAP,
    "input_parts": input_parts, "input_content_hashes_computed": False,
    "input_kind": "synthetic_test" if header.get("SyntheticTest", 0) else "user_catalog",
    "ngroups_total": expected_groups, "nsubhalos_total": expected_subhalos,
    "base_sample_count": len(clean),
    "catalog_sha256": hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
    "python": platform.python_version(),
    "packages": {name: version(name) for name in ("numpy", "pandas", "scipy", "h5py",
                                                  "matplotlib")},
    "units": {"Mstar": "Msun", "MBH": "Msun (sum in subhalo)",
              "M200c_parent": "Msun", "Mdot_BH": "Msun/yr",
              "SFR": "Msun/yr", "positions": "ckpc/h", "lambda_edd": "dimensionless"},
    "note": "Derived base sample; lower mass cuts require reprocessing original input.",
}
(OUTPUT_DIR / "catalog_metadata.json").write_text(
    json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
)
print("Izvezeni su katalog, metapodaci, tabele i grafikoni u outputs/.")

```
