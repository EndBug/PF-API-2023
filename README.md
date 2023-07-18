# PF-API-2023

> This project is in Italian.

## Struttura

Questo repo contiene la mia soluzione alla prova finale del corso di Algoritmi e Principi dell'Informatica, AA 2022-2023.

Il sorgente è in [`./solution.c`](./solution.c), e può essere compilato con `make solution` (oppure `just b`).  
Altri script utili per il testing si trovano in [`./justfile`](./justfile) (usando [`casey/just`](https://github.com/casey/just)).  
La traccia della prova, assieme ai tool consigliati, si trova in [`./docs/`](./docs/), mentre i test case e il generatore si trovano in [`./test/`](./test/).

## Idea

Le stazioni sono salvate sotto forma di albero rosso-nero, così da ottimizzare ricerca, inserimento, ed eliminazione.  
I veicoli di ciascuna stazione sono salvati in un max-heap, all'interno di ciascun nodo stazione. A parte aggiungere e rimuovere i veicoli, l'unica altra operazione che si svolge sui veicoli è la ricerca del massimo.  

Ogni volta che viene richiesto di calcolare un percorso, viene costruito un grafo contenente tutte le stazioni comprese tra la partenza e l'arrivo.  
A causa delle proprietà del problema, il grafo è diretto e aciclico ([Directed Acyclic Graph](https://en.wikipedia.org/wiki/Directed_acyclic_graph)), quindi si può utilizzare l'algoritmo di pathfinding per questo tipo di grafi per trovare il percorso più breve, con le seguenti accortezze:
- L'ordinamento topologico deve rispecchiare il verso di percorrenza della strada
- Il peso di tutti gli archi è unitario, dato che contiamo solo il numero di step
- Vengono considerati solo gli archi nel verso di percorrenza
- Per ottimizzare lo spazio, non vengono salvate né liste, né matrici di adiacenza, dato che si possono ottenere i nodi adiacenti a un altro scorrendo in ordine
- Si deve modificare la funzione di rilassamento per favorire sempre i nodi con chiave minore

Una volta terminato l'algoritmo, basta ripercorrere il percorso a ritroso partendo dalla destinazione.

Il grafo viene poi deallocato, mentre l'albero non viene alterato dal comando di ricerca del percorso.
