# Pregled projekta i plan istraživanja

Pregledan je `MilosP123/OAP2`, main commit
`64f1807ff55681d23a1e0188af11a2bb5023dc93` (17. septembar 2026).
Ovaj deo beleži stanje pre ispravki. Aktuelna revizija od 20. septembra 2026.
primenjuje tehničke ispravke i proširuje automatske provere, uz isti naučni fokus.
Originalni stvarni podaci nisu dostupni u radnom okruženju i nisu menjani.
Aktuelni status testova je u VERIFIKACIJA.md; konkretan kod pre/posle u IZMENE.md.

## Zatečena organizacija

| Fajl | Sadržaj |
|---|---|
| Projekat/Untitled7.ipynb | Glavna AGN analiza: 33 ćelije, sve tipa code, nekoliko praznih i zakomentarisanih |
| Projekat/README.md | Prazan (jedan newline) |
| Drugi_domaci_MP.ipynb | NASA Exoplanet Archive, vrući Jupiteri, bootstrap, GMM/AIC, ECDF |
| Treci_domaci_MP.ipynb | gr.dat, poređenje raspodela fluksa, KS/U test i simulacija ponašanja t-testa |
| Untitled3.ipynb | sunspot.npz, B i uglovi inc/az, statističke raspodele i fitovanje |
| gr.dat | Mala tabela histogramskih brojanja G/R |

U polaznoj grani nema requirements.txt, testova, AGENTS.md ni nastavnih .md
uputstava sa slika. `sunspot.npz` nije u repozitorijumu. Domaći zadaci nisu
menjani; detaljna revizija odnosi se na AGN projekat.

Glavni notebook radi samo sa jednim katalogom za z približno 0. Komentar o
z=0,1,2 nije implementacija višesnapshot analize. Učitava tri Group polja
(Group_M_Crit200, GroupFirstSub, GroupNsubs) i sedam Subhalo polja
(SubhaloFlag, SubhaloGrNr, SubhaloMassType, SubhaloBHMass, SubhaloBHMdot,
SubhaloSFR, SubhaloPos). Ne koristi pojedinačne BH čestice ili morfologiju.

## Sačuvani rezultati: samo istorijski izlazi

| Okruženje | Galaksije | AGN | Frakcija |
|---|---:|---:|---:|
| field | 14715 | 2610 | 0.177370 |
| group | 6029 | 320 | 0.053077 |
| cluster-like | 308 | 12 | 0.038961 |

Ukupno je prikazano 21052 galaksije i 2942 AGN kandidata. Korelacija
logM200c–loglambda ima sačuvani Pearson r=-0.491683 i Spearman rho=-0.392167,
na 18981 paru. Ovo nisu ponovo izračunati niti validirani rezultati.
Izlazi postoje čak i u ćelijama čiji je trenutni kod u celosti zakomentarisan.
Zato se ne mogu koristiti kao dokaz uspešnog Restart & Run All.

## Problemi utvrđeni u originalu (istorijski pregled)

1. KRITIČNO: ćelija 12 koristi nedefinisani M200c_all umesto M200c_a.
   NameError reprodukovan nezavisno, na malom nizu. Stari kernel može sakriti
   grešku i koristiti zaostalu verziju niza.
2. KRITIČNO: običan sorted sortira .10 pre .2. Globalni SubfindID i Group
   indeksi tada nisu očuvani za nepopunjene brojeve delova. Efekat na konkretne
   rezultate zahteva proveru stvarnih imena i ponovno računanje.
3. KRITIČNO: nema provere svih delova i potpunosti polja; ucitaj_polje tiho
   preskače nedostajuće polje. Proizvoljan podskup delova ne daje katalog čiji
   globalni indeksi počinju od nule i automatski ostaju važeći.
4. VISOKO: NaN lambda >= prag daje False; nedostajuće vrednosti postaju
   navodno neaktivne galaksije. Inf BH/Mdot i Mstar nisu potpuno kontrolisani.
5. VISOKO: centralne i satelitske galaksije nisu klasifikovane iako se
   GroupFirstSub učitava. GroupNsubs takođe nije korišćen za validaciju.
6. VISOKO: sirove AGN frakcije mešaju razlike u Mstar, statusu centrale/satelita
   i okruženju. Više galaksija iz istog halo-a nisu nezavisne realizacije okruženja.
7. VISOKO: apsolutna putanja i odsustvo podataka/preuzimanja blokiraju reprodukciju.
8. SREDNJE: ETA=0.1 je nedokumentovana postprocesna pretpostavka različita od
   TNG modelske efikasnosti 0.2; potrebno je ispitati osetljivost, ne tiho menjati.
9. SREDNJE: field je maseni razred, ne potvrđena izolovanost. Paul et al. daje
   motivaciju za približno 8e13 Msun; 5e12 je donja granica njihovog uzorka,
   ne univerzalna fizička granica izolovanih galaksija.
10. SREDNJE: nema neizvesnosti frakcija, izveštaja o selekciji po koracima,
    provere rezolucije ili dokumentovane definicije populacije.
11. NIŽE: ponovljeni grafikoni, neiskorišćeni importi, prazne ćelije, izostanak
    objašnjenja u Markdown ćelijama i sačuvani zastareli izlazi.

## Fizički izbori koji zahtevaju obrazloženje

Konverzije masa (puta 1e10/h) i akrecije (puta 1e10/0.978e9) već su ispravne.
Promena u float64 pre računanja je numeričko poboljšanje, ne nova fizička definicija.
SFR je već u Msun/yr. Koordinate treba označiti ckpc/h; fizička rastojanja
dobijaju faktor a/h i moraju uzeti u obzir periodične granice kutije.

SubhaloFlag==1 je opravdan osnovni rez. Mstar>=1e9 Msun je početni izbor,
ali sama vrednost ne dokazuje kompletnost ili konvergenciju AGN statistike.
Potrebno je proveriti rezultate za 1e9, 10**9.5 i 1e10 Msun, broj zvezdanih
čestica kada je dostupan, BH occupation fraction i uticaj modela formiranja BH.
Prag M200c>0 odbacuje neupotrebljivu masu; treba prijaviti broj takvih objekata.

Odnos dobijen iz zbirnih BH polja nije u opštem slučaju odnos jedne nuklearne
crne rupe: za jednako ETA predstavlja masom ponderisan prosek pojedinačnih
Eddingtonovih odnosa. Za izbor dominantnog nuklearnog BH potrebni su dodatni
podaci o BH česticama i jasan kriterijum izbora. Ne menjati postojeću definiciju
bez potrebe: za početak je tačno imenovati i navesti ograničenje.

Za lambda=0 ne računati logaritam, ali zadržati galaksiju u imeniocu AGN frakcije.
Posebno prikazati frakciju među svim klasifikovanim galaksijama i među onima
sa rezolviranom BH. Broj nepoznatih statusa prikazati po okruženju.

P-vrednost ispisana kao 0.0 nije dokaz apsolutne sigurnosti: numerička vrednost
može biti ispod opsega reprezentacije. Korelacije ne dokazuju uzročnost.
Korelacija lambda sa MBH takođe deli MBH kroz definiciju lambda proporcionalno
Mdot/MBH; potrebna je oprezna interpretacija i zbog matematičke zavisnosti.

## Predloženi grafikoni i statistika

| Korak | Grafik/statistika | Svrha i provera |
|---|---|---|
| Selekcija | Broj objekata posle svakog reza, razdvojeno po okruženju | Otkriva gde nastaje razlika uzoraka; brojevi moraju biti usklađeni |
| Populacije | Histogram/ECDF logMstar i udeo satelita po okruženju | Proverava da li se upoređuju slične galaksije |
| BH prisustvo | Udeo MBH>0 prema Mstar i okruženju | Razdvaja prisustvo BH od njene aktivnosti |
| Glavni rezultat | AGN frakcija prema Mstar, tri okruženja, odvojeno centrale/sateliti | Pokazuje ostaje li razlika pri kontroli mase i tipa |
| Kontinuirano okruženje | AGN frakcija prema logM200c sa brojevima galaksija i halo-a | Proverava da zaključak nije posledica proizvoljnih granica kategorija |
| Akrecija | ECDF pozitivne lambda i poseban udeo nulte akrecije | Izbegava da logaritmi prikriju neaktivne galaksije |
| BH/akrecija | 2D hist MBH–lambda i MBH–Mdot, iste ivice binova za poređenja | Razdvaja raspodele masa i akrecije; obeležiti prag AGN |
| Robusnost | AGN pragovi 1e-4,1e-3,1e-2; ETA=0.1/0.2; maseni rezovi | Proverava stabilnost znaka i veličine efekta |

Za frakcije prvo dati Wilson intervale kao binomnu referencu. Za glavno
zaključivanje koristiti bootstrap celih parent halo grupa (sa fiksnim seed-om),
uz kontrolu Mstar i statusa centrale/satelita. Prijaviti i broj jedinstvenih
halo-a; sa vrlo malo klastera ni bootstrap nije pouzdan. Prostorna korelacija
između halo-a i konačna zapremina ostaju dodatna ograničenja.

Fisher/hi-kvadrat za AGN–okruženje i KS/Mann–Whitney za raspodele mogu poslužiti
kao početne provere samo uz uvažavanje nezavisnosti. Ne permutovati proizvoljno
pojedinačne satelite između okruženja. Za poređenja po parovima korigovati
višestruka testiranja (npr. Holm). Uvek prikazati veličinu razlike i interval,
a ne samo p-vrednost. Spearman zadržati kao deskriptivnu dopunu.

Opcionalno: logistička regresija AGN statusa na logMstar, okruženje i tip
galaksije, sa halo bootstrap-om za neizvesnost. Ne unositi istovremeno masene
kategorije i istu halo masu bez jasnog razloga. Dodavanje MBH menja pitanje
u efekat pri fiksnoj BH masi; ono može ukloniti deo fizičke veze koju istražujemo.

## Aktuelna revizija i zahtevi profesorke

Primenjene su provere kompletnosti, dimenzija, tipova i indeksa; M200c_parent je
jedinstven naziv. Centrale i sateliti se određuju pre filtriranja. NaN/inf BH
status ostaje nepoznat, a nulta akrecija ostaje u imeniocu AGN frakcije.
Prikazana je tabela rezova, isti binovi za poređenje histograma, bezbedni prazni
grafikoni, frakcije po masi/tipu i halo bootstrap CI sa seed=42. Dodat je pregled
osetljivosti za ETA=0.1/0.2 i različite masene/AGN pragove. Opcionalni modeli
regresije i testovi raspodela nisu neophodni za tehničku ispravnost i ostaju
istraživačko proširenje, uz prethodnu proveru pretpostavki na stvarnim podacima.

| Zahtev | Aktuelno stanje |
|---|---|
| Relativne putanje, svi importi na početku | Ispunjeno u glavnom notebook-u |
| Bez interaktivnog unosa, instalacije i samoprovere u notebook-u | Ispunjeno; fizička validacija ulaza ostaje deo učitavanja |
| Fiksan seed | 42 za halo bootstrap |
| README i requirements | Dokumentacija i verzije direktnih zavisnosti dodate |
| Automatske provere | tests/test_catalog.py; fresh-kernel alat; GitHub Actions konfiguracija |
| Stari izlazi | Uklonjeni iz notebook-a; istorijske vrednosti iznad su jasno označene |
| Restart & Run All i pylint | Tačan status u VERIFIKACIJA.md |
| Dostupnost stvarnih podataka za drugi računar | Preostaje; veliki katalog je na korisnikovom disku |
| Dodatna nastavna .md uputstva | Nisu u repozitorijumu; provera je prema dostavljenim slikama |

## Preostali koraci koji zahtevaju stvarne podatke

1. Pokrenuti notebook nad svim delovima pravog TNG100-1 group kataloga.
2. Pregledati tabelu selekcije i brojeve nepoznatih vrednosti, ponovo izračunati
   AGN frakcije i uporediti ih sa istorijskim brojevima bez prilagođavanja koda
   da bi dao stare rezultate.
3. Dostaviti outputs/catalog_base.csv.gz i catalog_metadata.json. Ne slati samo AGN redove.
4. Proveriti stvarne raspodele Mstar, broj halo-a i stabilnost efekta u tabelama
   mass_fractions/sensitivity; bootstrap ne rešava malu kosmološku zapreminu.
5. Za ispit priložiti stvarni kompaktan ulaz i samostalni analitički notebook,
   ili pouzdano preuzimanje, pa izmeriti vreme ćelija pri timeout-u 60 s.
6. Proveriti ciljnu instalaciju i demonstrirati Restart & Run All uživo.

## Predložena struktura istraživačkog rada

1. Pitanje: razlikuje li se AGN aktivnost sa okruženjem pri uporedivoj Mstar
   i tipu galaksije? Razlikovati povezanost od uzročnog uticaja.
2. Podaci: TNG100-1, snapshot 99, polja, jedinice, poreklo i ograničenja rezolucije.
3. Metod: izbor uzorka i tabela selekcije, definicije BH/AGN/okruženja.
4. Rezultati: osnovna demografija, AGN frakcije pri kontroli mase, centrale/sateliti,
   raspodele akrecije i intervali neizvesnosti.
5. Robusnost: pragovi, ETA, BH occupation, broj halo-a i konačna zapremina.
6. Diskusija: mogući fizički mehanizmi uz ograničenja zbirnih BH polja,
   trenutnih akrecionih stopa i simulacionog feedback modela.
7. Zaključci: šta podaci podržavaju, veličina efekta, šta ostaje neodređeno.
8. Reprodukcija: notebook, ulazni izvedeni katalog, verzije i testovi.

Za 10-minutnu prezentaciju: 1 min pitanje/potreba, 1 min podaci, 2 min metod,
2 min ključni rezultat, 2 min živo pokretanje i objašnjenje ucitaj_polje ili
mapiranja indeksa, 1 min naučni izazovi, 1 min naučeno i ograničenja.
Sačuvani izlazi ne zamenjuju demonstraciju izvršavanja.

Izvori: [TNG dokumentacija](https://www.tng-project.org/data/docs/specifications/),
[Pillepich et al.](https://arxiv.org/html/1703.02970v2),
[Paul et al.](https://arxiv.org/pdf/1706.01916),
[SciPy pearsonr](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html).
