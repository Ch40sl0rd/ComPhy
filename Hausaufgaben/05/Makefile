# make
# ./hausaufgabe1
# Lars Döpper
#Compiler für Module und Programme
CC := gcc -g -ggdb -Wall -Wpedantic -std=c99
#Linker für fertige Programme
LD := gcc
#benötigte Bibliotheken
BIBS := -lm -lgsl -lblas

# alle Programmnamen, hier per Hand aufgelistet, muss manuell genändert werden
PROGRAMS := hausaufgabe5

# alle Quelldateien im momentanen Verzeichnis werden aufgelistest per "wildcard"
SOURCES := $(wildcard *.c)

# Alle zu kompilierenden Objektdatien, heir wird die Dateiendung
# von .c nach .o umgewandelt
OBJECTS := $(SOURCES:.c=.o)

# Alle Headerdateien
HEADERS := $(wildcard *.h)

# Die Programm- und Modul-objekte muessen unterschiedlich behandelt werden,
# wir filtern also alle Programmobjekte aus der Liste der Module
MODULES := $(filter-out $(addsuffix .o,$(PROGRAMS)), $(OBJECTS) )

# die erste Regel ist die Standardregel. Wir wollen, dass alle Programme
# kompiliert und gelinkt werden und nennen diese Regel konventionell "all"
all: $(PROGRAMS)

# ein praktisches Target, um die Variablen ausgeben zu lassen
list-variables:
	@echo HEADERS=$(HEADERS)
	@echo OBJECTS=$(OBJECTS)
	@echo MODULES=$(MODULES)
	@echo PROGRAMS=$(PROGRAMS)

# Regel zum Verlinken der ganzen Programme
$(PROGRAMS): % : %.o $(MODULES)
	$(LD) -o $@ $^ $(BIBS)

# die Programmobjekte haben andere Abhaengigkeiten als die Modulobjekte
# im letzten Feld machen wir aus programmN.o programmN.c um auf die Quelltextdatei
# fuer das Programmobjekt zu verweisen
# Die "SECONDEXPANSION" und die doppelten Dollarzeichen sind notwendig, 
# damit die Variablen und Funktionen
# auf der rechten Seite auch wirklich ausgewertet werden koennen...
# https://stackoverflow.com/questions/25611945/getting-the-target-name-in-a-function-in-prerequisites-list-makefile
.SECONDEXPANSION:
$(addsuffix .o, $(PROGRAMS)): % : $$(addsuffix .c,$$(basename $$(@))) $$(HEADERS) Makefile
	$(CC) -c -o $@ $< $(BIBS)

# Regel zum Kompilieren von Objekten, mit Ausnahme der Programmobjekte, dafuer
# haben wir ja eine spezielle Regel
%.o: %.c Makefile
	$(CC) -c -o $@ $< $(BIBS)

.PHONY: clean all list-variables
clean:
	rm $(OBJECTS) 

