Generatore di soluzioni pubblico Prova Finale API 2023
******************************************************

Per eseguire il generatore di soluzioni, seguire le istruzioni specifiche per
il vostro sistema operativo.

# Windows

Sono supportati Windows 10 e le versioni successive per processori Intel a
64 bit.

1. Aprire il prompt dei comandi
2. Digitare "cd" (senza virgolette) seguito da uno spazio e dal percorso dove
   si trova questo file. Per esempio, se questo file si trova nella cartella
   Download, e la vostra cartella utente è C:\Users\Studente,
   dovete digitare "cd C:\Users\Studente\Download".
3. Digitare "test_gen_2023_w64.exe" per eseguire il generatore. I comandi
   generati verranno stampati sul terminale. Per scriverli a file invece
   digitate "test_gen_2023_w64.exe > nomefile.txt", sostituendo "nomefile.txt"
   col nome file desiderato.
   
# Linux

È supportata Ubuntu 22.04 per processori Intel a 64 bit, il tool dovrebbe
funzionare anche su altre distribuzioni.

1. Aprire il Terminale della vostra distribuzione.
2. Digitare "cd" (senza virgolette) seguito da uno spazio e dal percorso dove
   si trova questo file. Per esempio, se questo file si trova nella cartella
   Downloads, e la vostra cartella utente è /home/studente,
   dovete digitare "cd /home/studente/Downloads".
3. Digitare "test_gen_2023_linux" per eseguire il generatore. I comandi
   generati verranno stampati sul terminale. Per scriverli a file invece
   digitate "test_gen_2023_linux > nomefile.txt", sostituendo "nomefile.txt"
   col nome file desiderato.
   
# macOS

È supportato macOS Big Sur 11.0 e le versioni successive per processori Intel
a 64 bit e per processori Apple Silicon.

1. Aprire Terminale.app
2. Digitare "cd" (senza virgolette) seguito da uno spazio e dal percorso dove
   si trova questo file. Per esempio, se questo file si trova nella cartella
   Downloads, e la vostra cartella utente è /Users/studente,
   dovete digitare "cd /Users/studente/Downloads".
4. Digitare "xattr -c test_gen_2023_macos". Questo comando elimina l'attributo
   quarantena aggiunto automaticamente ai file scaricati da internet che
   di default impedisce l'esecuzione di eseguibili non firmati da Apple
   (siccome siamo tirchi non intendiamo pagare 100$ a Apple per firmare un
   singolo eseguibile).
3. Digitare "test_gen_2023_macos" per eseguire il generatore. I comandi
   generati verranno stampati sul terminale. Per scriverli a file invece
   digitate "test_gen_2023_macos > nomefile.txt", sostituendo "nomefile.txt"
   col nome file desiderato.

# Opzioni disponibili

È possibile specificare opzioni aggiuntive al generatore aggiungendole
alla riga di comando. Le opzioni aggiuntive seguono il seguente formato.
Notare che ogni opzione può essere specificata con un nome lungo o corto.
Per esempio, per specificare il numero di stazioni è equivalente scrivere
"test_gen_2023_xxx -s 100" oppure "test_gen_2023_xxx --stations=100".

-k KM, --max-km=KM   Massimo valore chilometrico. Default: 100000000.
-s N, --stations=N   Numero medio di stazioni. Default: 30.
-c N, --min-cars=N   Numero minimo di auto per stazione. Default: 5.
-C N, --max-cars=N   Numero massimo di auto per stazione. Default: 10.
-r N, --routes=N     Numero medio di richieste di percorso, minimo 1.
                     Default: 10.
-d N, --removals=N   Rapporto percentuale tra comandi di rimozione e comandi
                     di aggiunta. Deve essere tra 0 e 49.
                     Default: 15.
-z N, --stupidity=N  Frazione percentuale di comandi errati. Deve essere tra
                     0 e 99. Default: 5.
