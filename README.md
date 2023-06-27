# PF-API-2023

> This project is in Italian.

## Idea

## Strutture dati

Le stazioni vengono salvate all'interno di un **RB tree**, dove le chiavi sono il km a cui si trova la stazione.  
I veicoli di ciascuna stazione sono a loro volta salvati in un **max-heap**, identificati dall'autonomia del veicolo.

```c
Station {
  int key;
  Station *left, *right;

  Vehicle *vehicles;
}

Vehicle {
  int range;
  Vehicle *left, *right;
}
```

## Algoritmo di ricerca

Per ogni stazione che considero guardo il veico con l'autonomia maggiore (per questo uso uno heap).
- Se l'autonomia è sufficiente per raggiungere la destinazione, ho trovato il percorso più breve.
- Altrimenti, per ogni stazione che posso raggiungere con l'autonomia massima ripeto ricorsivamente la stessa operazione. Poi confronto i risultati e prendo il migliore.
