# Okruženje galaksija i AGN aktivnost u TNG100-1

Projekat ispituje povezanost mase matičnog halo-a sa aktivnošću crnih rupa,
uz poređenje galaksija slične zvezdane mase i razdvajanje centralnih i satelitskih
galaksija. Glavni notebook je `Untitled7.ipynb`, snapshot 99 (z približno 0).
Analiza prečki/bar-a nije deo projekta.

**Stvarni ulazni podaci još nisu priloženi.** Notebook ne preuzima niti izmišlja
zamenske podatke. Rezultati su naučni tek nakon uspešnog izvršavanja nad stvarnim
TNG100-1 katalogom. Sintetički podaci u testovima proveravaju samo kod.

## Instalacija i pokretanje

Potrebni su Python 3.11 ili 3.12 i Jupyter. Iz direktorijuma `Projekat`:

```bash
python -m pip install -r requirements.txt
```

Paketi imaju zabeležene verzije. Instalaciju pokrenuti u istom okruženju koje
koristi Jupyter kernel. Potom otvoriti `Untitled7.ipynb` i izabrati
**Kernel → Restart & Run All**. Notebook ima sve importe u prvoj kodnoj ćeliji;
Markdown objašnjava postupak. Ne sadrži `input()`, `pip install`, `drive.mount`,
interaktivno raspakivanje ili ćelije za pokretanje testova.

## Stvarni podaci

Sve delove group kataloga TNG100-1, snapshot 99, smestiti u:

```text
data/group99/
```

Podržana su imena `groupcat-99.0.hdf5`, `groupcat-99.1.hdf5`, ... ili
`fof_subhalo_tab_099.0.hdf5`, ... . Ne koristiti obe kopije istog dela.
Numeracija mora pokriti 0 do `Header/NumFiles - 1`. Nedostajući deo ili polje
u nepraznoj tabeli zaustavlja analizu jasnom porukom; ne preskače se tiho.

Podaci mogu ostati na eksternom disku: lokalno napraviti simboličku vezu
`data/group99` ka direktorijumu stvarnih delova. Putanja u notebook-u ostaje
relativna. Veza ka privatnom disku nije rešenje za profesorkin računar i ne
predaje se na GitHub-u. Alternativa je kopiranje potrebnih fajlova u taj direktorijum.
Originalni HDF5 fajlovi se otvaraju isključivo u režimu `r`.

Notebook proverava h=0.6774, BoxSize=75000 ckpc/h i z≈0; to samo po sebi ne
potvrđuje nivo rezolucije. Korisnik mora obezbediti da je preuzeti izvor TNG100-1,
a ne TNG100-2/3. Manifest beleži nazive/veličine delova i njihove brojeve redova,
ali ne tvrdi da je izračunat hash kompletnog višegigabajtnog ulaza.

## Metodologija i jedinice

| Polje | Obrada/značenje |
|---|---|
| SubhaloMassType[:, 4] | Ukupna vezana zvezdana masa; množenje sa 1e10/h za Msun |
| SubhaloBHMass | Zbir masa BH po subhalou; množenje sa 1e10/h |
| SubhaloBHMdot | Zbir trenutnih akrecionih stopa; množenje sa 1e10/0.978e9 za Msun/yr; h se poništava |
| Group_M_Crit200 | Masa halo-a u sferi srednje gustine 200 puta kritična; množenje sa 1e10/h |
| SubhaloGrNr | Globalni indeks parent halo-a; ostaje isti nakon filtriranja |
| GroupFirstSub | Globalni indeks centralnog subhaloa; -1 za grupu bez subhaloa |
| SubhaloSFR | Već u Msun/yr |
| SubhaloPos | Izvorni ckpc/h; fizički kpc dobija se faktorom a/h |

Izvor definicija: [TNG specifikacije](https://www.tng-project.org/data/docs/specifications/).
Mstar nije aperturska masa. Zbirna BH polja ne identifikuju nužno jednu nuklearnu BH.

Osnovni uzorak `clean`: SubhaloFlag=1, konačna Mstar ≥ 1e9 Msun i konačna
M200c_parent > 0. Status centrale/satelita određuje se pre selekcije i ne menja
se ako centrala kasnije ne prođe rez. Tabela `selection` broji svaki korak.

Eddingtonov odnos računa se iz Lbol=ETA*Mdot*c² i Ledd=1.26e38*(MBH/Msun) erg/s.
ETA=0.1 je zadržana postprocesna pretpostavka. Nije izmerena efikasnost niti
univerzalno fizičko preslikavanje pri slaboj akreciji; TNG modelska vrednost je
0.2 ([opis modela](https://arxiv.org/html/1703.02970v2)). Tabela osetljivosti proverava oba izbora.

AGN kandidat: lambda_edd ≥ 1e-3. MBH>0 i Mdot=0 daje validnu lambda=0 i neaktivni
status. MBH=0 i Mdot=0 predstavlja galaksiju bez AGN-a u simulacionom katalogu.
NaN/inf, negativna akrecija ili nekonzistentni BH podaci daju nepoznat status.
Neaktivne galaksije se ne izbacuju iz imenioca. Logaritmi obuhvataju samo pozitivne
konačne vrednosti, pa njihove medijane nisu medijane svih galaksija.

Okruženja: field za M200c<5e12; group za 5e12≤M200c<8e13; cluster-like iznad toga,
u Msun. `field` je maseni razred, ne dokaz izolovanosti. Motivacija za gornji
prelaz je [Paul et al. 2017](https://arxiv.org/pdf/1706.01916); donja granica njihovog
uzorka nije univerzalna granica izolovanih galaksija.

## Statistika i izlazi

Glavna frakcija je N_AGN/N_poznatih. Tabele daju i N_ukupno, N_nepoznatih,
N_halo i moguće granice frakcije za celu populaciju. Te granice nisu intervali
poverenja. Dodatno se prikazuje frakcija među galaksijama sa MBH>0.

Maseni binovi imaju širinu 0.5 dex. Frakcije su odvojene po okruženju i
centrala/satelit statusu. Percentilni 95% intervali dobijaju se sa 1000 bootstrap
realizacija celih halo-a, seed=42. Za <2 halo-a interval nije definisan;
za mali broj halo-a i degenerisane uzorke može biti nepouzdan. Bootstrap ne
obuhvata kosmičku varijansu i korelacije između različitih halo-a. Poređenje
binova nije potpuno izjednačavanje masenih raspodela unutar svakog bina.

Pearson/Spearman su deskriptivni. Njihove naivne p-vrednosti nisu kontrolisani
test uticaja okruženja niti dokaz uzročnosti. P=0 je numerički ispis.

Posle uspešnog izvršavanja, `outputs/` sadrži:

| Fajl | Namena |
|---|---|
| catalog_base.csv.gz | Osnovni uzorak, sa neaktivnim i nepoznatim statusima |
| catalog_metadata.json | Pragovi, jedinice, izvor, seed, verzije, SHA-256 izvedenog kataloga |
| selection.csv | Broj objekata posle svakog reza, po okruženju |
| summary.csv / summary_central.csv / summary_bh.csv | Demografija, imenilac i AGN frakcije |
| correlations.csv | Korelacije i broj konačnih parova |
| mass_fractions.csv | AGN frakcije i halo bootstrap intervali po masi/tipu/okruženju |
| sensitivity.csv | ETA 0.1/0.2, AGN pragovi 1e-4/1e-3/1e-2, maseni rezovi 1e9/10^9.5/1e10 |
| lambda_edd_vs_m200c.png | Pozitivna akrecija prema masi halo-a i AGN prag |
| bh_raw.png / bh_clean.png / bh_clean_lognorm.png | 2D raspodele sa zajedničkim binovima |
| agn_fraction_vs_mstar.png | Glavno poređenje okruženja pri kontroli mase/tipa |

Za nastavak analize dovoljno je poslati `catalog_base.csv.gz` i
`catalog_metadata.json`. Ovaj izvoz ne omogućava spuštanje Mstar reza ispod
onog koji je korišćen pri izvozu. Ne slati samo AGN redove.

## Automatske provere

Iz `Projekat`:

```bash
python -m pytest -q tests
python tools/check_notebook.py
python -m nbqa pylint Untitled7.ipynb
```

Druga komanda kopira neizmenjeni notebook u privremeni direktorijum, generiše
sintetički katalog i pokreće pravi `pytest --nbmake --nbmake-timeout=60` u novom
kernelu. Privremeni ulazi/izlazi ne mogu prepisati stvarne podatke.
Testovi proveravaju i bajt-po-bajt nepromenjenost ulaznih HDF5 fajlova.

Nad stvarnim ulazom profesorkina provera glasi:

```bash
pytest --nbmake --nbmake-timeout=60 Untitled7.ipynb
nbqa pylint Untitled7.ipynb
```

GitHub Actions konfiguracija proverava Python 3.11/3.12, regresione testove,
pylint bez E poruka i ocenu najmanje 5, zatim fresh-kernel nbmake nad sintetičkim
podacima. To ne potvrđuje naučne rezultate niti brzinu čitanja velikog stvarnog
kataloga. Detaljni poslednji rezultati su u `VERIFIKACIJA.md`.

  objašnjenje funkcije, naučni izazovi, naučeno i ograničenja; 15 minuta pitanja.
- Na ciljnom računaru proveriti `git status` i da lokalni commit odgovara
  predatom GitHub commit-u. Čist working tree sam po sebi ne dokazuje da je urađen push.

