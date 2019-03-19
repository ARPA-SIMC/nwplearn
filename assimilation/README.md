## Optimal interpolation

### opt_int

Calcola l'analisi grigliata a partire da un background costante (=10.)
e una serie di osservazioni fornite dall'utente, ad esempio nei file
`obs_input_2.naml` `obs_input_3.naml`. Il file di osservazioni deve
essere fornito come parametro a linea di comando. Le caratteristiche
della griglia del background e dell'analisi sono lette dal file
`grid.naml`.

#### Esercizio

Verificare come cambia il risultato dell'analisi al cambiare dei
parametri delle matrici di varianza/covarianza dell'errore di
background e osservazioni.

Modificare la varianza dell'errore del background in `grid.naml` e la
varianza dell'errore delle osservazioni nel file di osservazioni
scelto.

Modificare scala spaziale della funzione analitica che definisce le
matrici di covarianza dell'errore del background e delle osservazioni
modificando le righe:

```
an_opt_int%o_mat_func = 'gauss'
an_opt_int%o_mat_fact = 0.1_fp_d
an_opt_int%b_mat_func = 'gauss'
an_opt_int%b_mat_fact = 0.4_fp_d
```

in `opt_int.F90` (questa seconda modifica richiede ricompilazione con
`make`).

Analizzare il risultato in `opt_int.pdf`.

### opt_int_multi

Calcola l'analisi grigliata a partire da un background costante (=10.)
e una serie di osservazioni via via più fitte che definiscono un campo
predefinito con funzioni di Fourier bidimensionali più opportune
perturbazioni casuali proporzionate all'errore osservativo fornito.

Le caratteristiche della griglia del background e dell'analisi sono
lette dal file `grid.naml`. Il risultato dovrebbe essere un campo di
analisi il più possibile vicino alla verità, pari alla funzione di
Fourier bidimensionale.

#### Esercizio

Ottenere una minimizzazzione dell'errore dell'analisi per diversi
numeri di osservazioni intervenendo sulle caratteristiche delle
matrici di varianza/covarianza dell'errore di background e
osservazioni

Modificare la varianza dell'errore del background in `grid.naml` e la
varianza dell'errore delle osservazioni nella riga di programma:

```
! observation error distribution variance
obs_set%err = 0.1_fp_d
```

Modificare la scala spaziale della funzione analitica che definisce le
matrici di covarianza dell'errore del background e delle osservazioni
modificando le righe:

```
an_opt_int%o_mat_func = 'gauss'
an_opt_int%o_mat_fact = 0.1_fp_d
an_opt_int%b_mat_func = 'gauss'
an_opt_int%b_mat_fact = 0.4_fp_d
```

in `opt_int_multi.F90` (queste ultime modifiche richiedono ricompilazione
con `make`).

Analizzare il risultato in `opt_int_multi.pdf`.

