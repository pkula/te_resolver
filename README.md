# te_resolver

## Goal:
Create pipeline to find TE within genome

## Dictionary:
TE - transposable elements
ONT - Oxford Nanopore Technologies

## Data:
Reference genome
Transcposon sequences (TE)
Long reads from ONT sequencing

## Steps
* Identification of ONT reads containing TE sequences
- blast TE do odczytów (-dust no -perc_identity 0.9) - tu można dodać dodatkowy etap filtrowania hitów o minimalnej długości równej np. połowie długości TE?
-  na podstawie pliku wynikowego tworzony jest plik z listą nazw, na podstawie którego wybierane są sekwencje ONT do nowego pliku i plik Bed z lokalizacją miejsca dopasowania do transpozonu, który jest wykorzystywany do maskowania sekwencji TE w odczycie. Plik Bed można, z koordynatami łączonymi -d = 1500 dla LTR-RT i 200 dla MITE) można przefiltrować po długości

* Identification of insertion sites in the reference genome
-blast odczytów z zamaskowanym rejonem do genomu referencyjnego (-perc_identity 0.8 -max_target_seqs 4 -out TE_sample_genome.bl)

