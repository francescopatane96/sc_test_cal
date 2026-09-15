# Bozza di manoscritto e figure

Titolo di lavoro neutro: **Cellular composition and immune transcriptional states associated with established PIRA in cerebrospinal fluid**. Specificare la coorte femminile nel titolo finale o nell'abstract. Adattare il titolo ai risultati reali; nessuna firma o meccanismo è già dimostrato.

## Domanda centrale

Nei soggetti con PIRA già avvenuta, il CSF presenta differenze nella composizione immunitaria e negli stati trascrizionali, e queste associazioni persistono considerando età ed EDSS alla diagnosi? L'informazione trascrittomica migliora la classificazione retrospettiva oltre alle covariate cliniche?

| Figura | Pannelli suggeriti dalla pipeline | Messaggio da valutare, non risultato assunto |
|---|---|---|
| 1. Coorte e atlante | Diagramma esclusioni; tabella soggetti; distribuzioni Age/EDSS; UMAP L2 e mieloidi | Disegno, compartimenti, numerosità effettiva; PIRA precedente al campionamento |
| 2. Composizione | Proporzioni per soggetto; effetti aggiustati; confronto L1/L2 e mieloidi | Differenze riproducibili fra individui, distinguendole dai soli conteggi di cellule |
| 3. Stati e pathway | Mediane UCell, risultati CAMERA; pannelli focali interferoni, TNF/NFκB, complemento, metabolismo, presentazione antigenica quando supportati dai set | Programmi coordinati anche se i singoli DEG sono pochi |
| 4. Approfondimento mieloide | Pseudobulk DEG, ORA up/down, score manuali descrittivi, LOPO | Eventuale stato mieloide associato a PIRA, supportato da più soggetti |
| 5. Covariate e classificazione | RF clinica/genica/combinata su medesimi split; incremento AUC; score paziente; SHAP fuori training | Valore aggiunto rispetto ad Age ed EDSS, oppure sua assenza |
| 6. Validazione, se disponibile | Coorte indipendente; citometria/proteine/validazione funzionale mirata; replica per soggetto | Conferma del risultato prioritario e sua interpretazione biologica |
| Supplementari | QC, missing, composizione sangue/CSF controlli, singola cellula vs pseudobulk, SCPA, modelli non aggiustati, stabilità e risultati negativi | Trasparenza, dipendenza da individui e limiti dell'inferenza |

Questa è una scaletta editoriale, non un impegno a sei figure principali. Se i risultati più solidi sono due, concentrare il manoscritto su quelli. La sola moltiplicazione delle analisi non rende uno studio adatto a Nature Communications: servono una domanda convincente, robustezza e una validazione proporzionata alla conclusione. Le figure di validazione non sono generate dalla pipeline in assenza di dati.

## Abstract da completare, senza risultati inventati

**Background.** La progressione indipendente dalle ricadute può associarsi a stati immunitari del compartimento liquorale non descritti dalla sola disabilità clinica.

**Methods.** Analisi esplorativa di scRNA-seq del CSF in donne con PIRA preesistente, RRMS e controlli, con un confronto indipendente CSF–blood nei controlli. Annotazioni Azimuth e mieloidi manuali; proporzioni e trascrittomica analizzate a livello di soggetto; modelli aggiustati per età ed EDSS alla diagnosi nei confronti fra pazienti. Classificazione retrospettiva con validazione a soggetti disgiunti.

**Results.** Inserire n soggetti e cellule dopo esclusioni, effetti con incertezza, FDR, stabilità per individuo e incremento OOF rispetto al modello clinico. Riportare anche assenza di incremento, se osservata.

**Conclusions.** Limitare le conclusioni ad associazioni nella coorte studiata. Non usare “predittore di futura PIRA” o “meccanismo causale” sulla base di questo disegno.

## Informazioni ancora da inserire nei Methods

Definizione operativa di PIRA e durata della conferma; date di diagnosi, EDSS, evento e prelievo; diagnosi dei controlli; ricadute e steroidi vicino al prelievo; trattamento; durata di malattia; attività RM; piattaforma, librerie e batch; QC/doublet/ambient RNA; origine dei counts; reference Azimuth e gene set delle annotazioni manuali. Questi dati non sono ricavabili dai nomi delle colonne forniti.

Descrivere le esclusioni per sesso come restrizione della popolazione analizzata e discutere la generalizzabilità. Non eliminare ulteriori pazienti in base a SHAP o alla capacità del modello di classificarli.

## Analisi ulteriori condizionate ai dati

Milo o altri test su neighborhood possono essere utili per stati continui se vi sono repliche sufficienti e una struttura di vicinato verificata; non sono implementati qui. Interazioni ligando–recettore sono ipotesi di comunicazione, da verificare senza sostituire i replicati biologici. VDJ richiede dati recettoriali; sopravvivenza richiede tempi e follow-up; traiettorie richiedono un processo e una radice plausibili. Non aggiungerli soltanto per ottenere un risultato significativo.

## Riferimenti metodologici e biologici di partenza

- [Schafflick et al., atlante sangue e CSF nella sclerosi multipla](https://www.nature.com/articles/s41467-019-14118-w): contesto per distinguere composizione e stati nei compartimenti.
- [muscat: confronti di stato in scRNA-seq multisoggetto](https://pmc.ncbi.nlm.nih.gov/articles/PMC7705760/): importanza della replicazione biologica.
- [scProportionTest](https://github.com/rpolicastro/scProportionTest): metodo richiesto, mantenuto come supporto cellulare esplorativo.
- [SCPA](https://jackbibby1.github.io/SCPA/reference/compare_pathways.html): confronto di distribuzioni di pathway senza selezionare una lista DEG.
- [UCell](https://www.bioconductor.org/packages/release/bioc/html/UCell.html): firme basate sui ranghi.
- [edgeR](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html) e [limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html): pseudobulk e gene set a livello di soggetto.

Il file `methods.md` di ogni modulo registra il contrasto e i parametri effettivamente usati; integrare la versione dei pacchetti da `sessionInfo.txt`. Conservare tutte le tabelle per rendere tracciabile la selezione dei pannelli finali.
