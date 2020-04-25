# make versucht standardmaessig das erste "Target" zu erstellen

# linking
programm: modul1.o modul2.o programm.o
	gcc -o programm modul1.o modul2.o programm.o
      
# kompilieren
modul1.o: modul1.c
	gcc -std=c99 -Wall -Wpedantic -c -o modul1.o modul1.c
      
modul2.o: modul2.c
	gcc -std=c99 -Wall -Wpedantic -c -o modul2.o modul2.c

programm.o: programm.c modul1.h modul2.h
	gcc -std=c99 -Wall -Wpedantic -c -o programm.o programm.c

