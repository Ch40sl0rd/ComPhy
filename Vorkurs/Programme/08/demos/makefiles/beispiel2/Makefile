# Kommandozeilenargumetne fuer gcc als Compiler
CC := gcc -std=c99 -Wall -Wpedantic
# als Linker benoetigt gcc nicht unbedingt Kommandozeilenargumente
LD := gcc

all: programm1 programm2

# beim Kompilieren '-c' nicht vergessen
programm1.o: programm1.c modul1.h modul2.h
	${CC} -c $<

programm2.o: programm2.c modul2.h modul3.h
	${CC} -c $<        

%.o : %.c %.h
	${CC} -c $<

# beim Linken '-o' nicht vergessen
programm1: programm1.o modul1.o modul2.o
	${LD} -o $@ $^

programm2: programm2.o modul2.o modul3.o
	${LD} -o $@ $^

