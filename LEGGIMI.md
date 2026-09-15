# Pipeline PIRA–RRMS, versione 3: esclusione dei cinque soggetti maschi

Questa versione usa **scc** come oggetto completo e **myeloid** per trasferire le annotazioni manuali, tramite barcode identici. `sample` identifica il soggetto. Non modifica gli oggetti originali in memoria: crea un oggetto di lavoro filtrato.

## Avvio in RStudio

1. Estrarre tutta la cartella dello ZIP, mantenendo insieme i file.
2. Caricare `scc` e `myeloid` nella stessa sessione R. `myeloid` deve contenere `myeloid_cluster` e `sample`; le sue cellule devono provenire da `scc`.
3. Se mancano pacchetti, aprire `INSTALL.R` ed eseguirlo. Installazione e download MSigDB richiedono internet. Non mescolare questi script con quelli precedenti.
4. Aprire `RUN_ALL.R` e premere **Source**. La cartella dei risultati è `PIRA_results_female_v3`, nella directory di lavoro R (`getwd()`).
5. Controllare `MODULE_STATUS.csv`, `MODULE_OUTPUT_STATUS.csv`, `subject_exclusions.csv`, `patient_metadata_female.csv` e, in ciascun confronto, `cluster_eligibility.csv`. Un modulo FAILED o non eleggibile NON è un risultato biologico negativo.
6. Per la Random Forest, aprire `RUN_RF.R` e premere Source nella stessa sessione, dopo il core. È separata perché costosa. Non basta eseguire singole righe della RF senza gli oggetti creati dal core.

Tutti i parametri sono in `00_config.R`; gli altri file sono caricati automaticamente. Nessun file del dataset è sovrascritto. Riutilizzare la stessa directory di output sovrascrive le esportazioni omonime: cambiare `cfg$output` per confrontare configurazioni.

## Scelte definite prima dei risultati

- Esclusione globale di **3611, 3553, 1743, 3608, 3627**, anche da controlli/blood se presenti. Identificativi numerici o factor sono convertiti in testo e ripuliti dagli spazi. La classificazione come femmine degli altri soggetti si basa sulla lista fornita, non su una colonna sesso verificata.
- Non viene aggiunto il sesso al modello dopo l'esclusione: sarebbe costante. Questa analisi riguarda quindi la coorte femminile. Se tutti i maschi sono PIRA, una regressione non crea l'informazione mancante sui maschi RRMS.
- PIRA è già avvenuta al prelievo: si studiano **associazioni con uno stato preesistente e classificazione retrospettiva**. EDSS è quello alla diagnosi, come indicato. Non si stima il rischio futuro di PIRA.
- PIRA vs RRMS: solo CSF, modello `~ Age + EDSS + group`. edgeR include anche `~ group` sugli stessi complete cases. Il confronto fra i due modelli aiuta a distinguere associazione totale e associazione condizionata alle covariate; non identifica causalità.
- PIRA vs ctrl e RRMS vs ctrl: solo CSF, modello `~ Age + group`. Non si inventa un EDSS uguale a zero nei controlli.
- Ctrl CSF vs ctrl blood: soggetti indipendenti, `~ Age + group`. Il codice si ferma su questo assunto se uno stesso ID controllo compare nei due tessuti, chiedendo di chiarire il pairing.
- L1, L2 e cluster mieloidi sono risoluzioni dello stesso dataset, non repliche indipendenti. I denominatori delle proporzioni sono tutte le cellule della rispettiva vista: le proporzioni mieloidi sono relative alle cellule di `myeloid`.
- Non viene stimata un'interazione malattia × tessuto in assenza di PIRA/RRMS blood. Il contrasto nei controlli descrive differenze di compartimento, non dimostra quali cambiamenti PIRA siano specifici del CSF.

## Tabella delle analisi implementate

Ogni modulo esportato contiene `results.csv`, `figure.png`, `figure.pdf` e `methods.md`. Un pannello vuoto documentato segnala assenza di test eleggibili. Alcuni moduli aggiungono dati sorgente e RDS. I Methods sono bozze da verificare e integrare con i dettagli clinici e sperimentali.

| Analisi | Unità e metodo | Parametri principali | Figura/output e interpretazione |
|---|---|---|---|
| Audit e atlante | Soggetti, cellule; UMAP esistente | Esclusione cinque ID; controlli metadata | Numerosità, UMAP L1/L2/mieloidi; audit esclusioni e counts |
| Proporzioni principali | Una proporzione per soggetto/pop./tessuto; limma asin-sqrt | Age; EDSS solo PIRA/RRMS; BH popolazioni | Boxplot con punti soggetto; effetto aggiustato e differenza descrittiva in punti percentuali |
| scProportionTest | Permutazioni e bootstrap sulle cellule aggregate | 10.000 repliche | Plot nativo; supporto esplorativo, non test indipendente sui pazienti |
| DEG pseudobulk | Somma RNA counts per soggetto/pop./tessuto; edgeR QL robust | ≥10 cellule/pop./soggetto; ≥4 soggetti/gruppo; TMM, filterByExpr | Volcano, tabelle adjusted/unadjusted; BH geni × cluster per contrasto/vista/modello |
| DESeq2 opzionale | Conteggi pseudobulk autentici interi | Stesso disegno; filtro ≥10 counts in ≥4 soggetti; poscounts | Volcano per cluster; sensibilità, non scelta del metodo più significativo |
| DEG singola cellula | Seurat Wilcoxon, log-normalizzazione RNA | min.pct 0,1; massimo 200 cellule/soggetto; nessun filtro FC iniziale | Volcano esplorativo; non risolve la pseudoreplicazione |
| Firme immunologiche | UCell per cellula, mediana per soggetto/cluster; limma | Hallmark; maxRank 1.500; ≥15 geni mappati | Score sorgente, confronto aggiustato e dotplot; score non equivale a funzione biologica misurata |
| Firme manuali esistenti | Mediane dei metadata `score_*` | Nessun test confermativo automatico | Grafici descrittivi; la firma usata per annotare un cluster non è validazione indipendente |
| Pathway senza soglia DEG | CAMERA su voom pseudobulk | Hallmark/Reactome, 15–500 geni; correlazione intergenica 0,01 | Dotplot direzionale; BH pathway × cluster per collezione |
| ORA up/down | Test ipergeometrico separato per direzione | BH DEG globale <0,05; abs(log2FC) ≥0,25; almeno 10 DEG | Dotplot, geni di overlap; universo = geni testati nel cluster; BH anche sulle due direzioni |
| SCPA di supporto | Confronto distribuzioni cellulari Hallmark | Se <10 DEG; circa ≤300 cellule/condizione, bilanciate per soggetto; 3 sottocampionamenti | Qval nativo, alto = maggiore perturbazione; non è Storey q-value e non aggiusta Age/EDSS |
| Correlazioni cliniche | Spearman entro condizione su mediana soggetto/cluster | Age/EDSS; ≥6 soggetti e ≥3 valori distinti | Dotplot rho, BH; esplorativo, senza ulteriore aggiustamento |
| PCA | LogCPM pseudobulk | Per popolazione e contrasto | Soggetti etichettati; utile per outlier e separazioni tecniche |
| Influenza LOPO | Ricalcolo edgeR omettendo un soggetto | Top 20 geni del modello adjusted | Dispersione log2FC; evidenzia dipendenza da individui, non validazione esterna |
| RF ripetuta | Un record per soggetto CSF; caret/ranger | 10 ripetizioni, fino a 5 fold stratificati; 300 alberi; min.node.size 3 | AUC, average precision, Brier, balanced accuracy; score OOF per soggetto |
| Confronto modelli RF | EDSS, Age, Age+EDSS, geni, firme, composizione, geni+covariate, completo | Identici split; 100 geni variabili selezionati solo nel training; nessun tuning sull'OOF | AUC appaiate e incremento rispetto ad Age+EDSS; nessun t-test tra ripetizioni dipendenti |
| SHAP | Solo soggetti tenuti fuori dal training; background training | Fastshap, 30 simulazioni; modelli geni/geni+covariate/completo | Importanza, copertura della selezione, plot direzionale; nessuna interpretazione causale |
| Permutazioni RF opzionali | Permutazione etichette soggetto, intera CV rifatta | 199 permutazioni; molto costoso | Distribuzione nulla e p corretto BH fra modelli; non è test di informazione genica condizionata ad Age/EDSS |

SCPA usa il formato di pathway documentato dagli autori e il suo Qval nativo: [documentazione SCPA](https://jackbibby1.github.io/SCPA/reference/compare_pathways.html). scProportionTest permuta le etichette delle cellule; per questo è accompagnato da un modello per soggetto: [codice del metodo](https://github.com/rpolicastro/scProportionTest/blob/master/R/permutation_test.R).

## Interpretazione dei dati e soglie

La matrice deve essere `RNA/counts`, con conteggi autentici. Il codice unisce i layer Seurat v5 e usa lo slot counts in v4. Non usa integrated, scale.data o esponenziazioni di dati corretti per ricostruire conteggi. La posizione nel layer counts da sola non dimostra la provenienza: verificarla se l'oggetto è stato modificato in precedenza. Conteggi frazionari vengono segnalati: edgeR li accetta, ma sono appropriati solo se rappresentano effettivamente conteggi stimati/corretti; DESeq2 viene saltato senza arrotondamenti arbitrari.

La FDR primaria è BH, non viene cambiata per ottenere positività. ORA può non produrre risultati; CAMERA è eseguito anche senza DEG. Il numero minimo di cellule non è una soglia di potenza: dopo le esclusioni molti cluster potrebbero avere pochi soggetti. Missing clinici e disegni non identificabili sono segnalati, senza riempire EDSS arbitrariamente.

La RF sui geni somma tutte le cellule CSF per soggetto: può catturare anche cambiamenti di composizione. I confronti pseudobulk per cluster aiutano a distinguere composizione e stato cellulare. Le probabilità RF non sono ancora calibrate per uso clinico. Le ripetizioni non aumentano il numero dei pazienti; i range fra ripetizioni non sono intervalli di confidenza. SHAP di un classificatore debole resta esplorativo. Feature correlate, incluso HLA, possono condividere o scambiarsi importanza.

Non sono aggiustati automaticamente trattamento, batch, durata di malattia o attività RM: non sono stati forniti. `seq_folder` non viene assunto essere batch indipendente dal soggetto. I geni mitocondriali/ribosomiali o HLA dominanti richiedono controlli di qualità e sensibilità prima di conclusioni biologiche.

Per aggiungere firme pubblicate, salvare un RDS contenente una lista nominata di vettori di simboli HGNC e indicare il percorso in `cfg$custom_signatures_rds`. Conservare per ogni firma riferimento, specie, direzione e versione; non costruirla dai DEG dello stesso confronto per poi presentarla come conferma indipendente.

## Verifica eseguita su questa versione

Parsing degli script R e test sintetici su esclusione degli ID (factor/numerici/spazi), disegno con Age/EDSS, assenza di EDSS nei controlli, collinearità, somme pseudobulk e universo ORA. **Non è stata eseguita la pipeline sul tuo `scc`**: l'oggetto non è disponibile in questo ambiente e mancano alcuni pacchetti analitici. L'esecuzione completa, i risultati biologici e la compatibilità delle versioni installate vanno verificati nella tua sessione R.
