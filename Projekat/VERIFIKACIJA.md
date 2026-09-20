# Verifikacija — 20. septembar 2026.

## Opseg i ulaz

Pregledana je aktuelna GitHub grana main na commit-u
`64f1807ff55681d23a1e0188af11a2bb5023dc93`, raniji PREGLED.md, IZMENE.md,
VERIFIKACIJA.md i kompletan pripremljeni notebook. Stvarni TNG podaci nisu
bili dostupni. Nisu menjani originalni podaci, domaći zadaci niti naučni fokus.

## Lokalno izvršene provere

| Provera | Rezultat |
|---|---|
| `MPLBACKEND=Agg python -m pytest -q tests` | 36 passed in 13.43s |
| Sve kodne ćelije redom u potpuno novom Python procesu | Prošlo, deo testa test_fresh_python_process |
| Notebook, standardni `nbqa pylint` | 9.42/10; bez E/F poruka, postoje C/R/W poruke |
| `python -m nbqa pylint Untitled7.ipynb --errors-only` | Exit 0 |
| `python -m nbqa pylint Untitled7.ipynb --fail-under=5` | Exit 0 |
| `python -m pylint tools/check_notebook.py tests/test_catalog.py --errors-only` | Exit 0 |
| JSON/schema i Python sintaksa sva četiri notebook-a | Prošlo; jedan SyntaxWarning u odvojenom domaćem zadatku |
| Pravi fresh-kernel nbmake nad sintetičkim ulazom, lokalno | Blokiran pre prve ćelije: mrežni soketi kernela nisu dozvoljeni (`Kernel died before replying to kernel_info`, `Operation not permitted`) |
| Stari izlazi u glavnom notebook-u | Uklonjeni, execution_count=None |
| Izvorni HDF5 fajlovi u testu | SHA-256 pre/posle isti |

Testovi koriste 12 sintetičkih delova sa 24 subhaloa, plus namerno neispravne
varijante i slučajeve bez akrecije/praznog uzorka. Header nosi SyntheticTest=1,
a izvoz metadata.input_kind=synthetic_test. Nijedan takav izlaz nije naučni rezultat.
Ne dozvoljava se automatska zamena nedostajućih stvarnih podataka sintetičkim.

Standardni pylint vraća exit 28 zbog C/R/W poruka pri podrazumevanom pragu 10;
posebno su potvrđeni uslov najmanje 5 i odsustvo E poruka. Nisu isključivane
E provere da bi se postigla ocena.

## Šta je provereno testovima

Numerički redosled delova; nedostajući/duplirani delovi; nedostajuća polja i
header atributi; pune dimenzije skalarnih i vektorskih polja; celobrojni tipovi;
negativni i preveliki parent indeksi; satelit pogrešno označen kao centrala;
prazni delovi; grupe bez subhaloa i signed/unsigned sentinel -1; neispravna
zaglavlja; konverzije masa/akrecije; nezavisna numerička provera lambda;
jednakost na AGN/okruženje pragovima; nulta i nepoznata akrecija; NaN/inf halo
mase; imenilac i granice AGN frakcije; prazan LogNorm grafikon; maseni maksimum
na granici bina; reproduktivan bootstrap halo-a; izvoz, kontrolna suma i
nepromenjenost ulaznih fajlova; novi Python proces bez prethodnih promenljivih.

## GitHub Actions

Dodata je konfiguracija `.github/workflows/agn-checks.yml` za Python 3.11/3.12.
Ona pokreće testove, oba pylint uslova i pravi `nbmake` preko
`tools/check_notebook.py`. GitHub rezultat se potvrđuje zasebno nakon objave;
sama konfiguracija nije dokaz da je provera izvršena ili prošla.

## Ograničenja koja ostaju

- Nema potvrde pravih AGN frakcija, broja galaksija ili fizičkih zaključaka
  bez stvarnog TNG kataloga. Istorijski brojevi su označeni u PREGLED.md.
- Nije izmereno trajanje učitavanja punog kataloga sa eksternog diska;
  ne može se potvrditi timeout=60 s za stvarne podatke.
- Podaci/preuzimanje za računar profesorke još nisu obezbeđeni.
- Originalna dodatna nastavna .md uputstva nisu u repozitorijumu; provereni
  su zahtevi vidljivi na četiri dostavljene slike.
- Domaći zadaci van Projekat/ su pregledani strukturno/sintaksno, ali nisu
  izvršeni niti menjani: sunspot.npz nedostaje, egzoplanetni notebook ima
  spoljašnji servis i zasebne zavisnosti. AGN requirements i CI važe za Projekat/.
- CI na Linux-u ne potvrđuje instalaciju na korisnikovom Mac-u.

## Testirano lokalno okruženje

Python 3.12.14

- numpy==2.3.5
- pandas==2.2.3
- scipy==1.17.0
- matplotlib==3.10.8
- h5py==3.16.0
- pytest==9.1.1
- nbmake==1.5.5
- nbqa==1.9.1
- pylint==4.0.8
- ipykernel==7.3.0
