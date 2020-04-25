# Compiler
CC=gcc -g -std=c99 -Wall -Wpedantic

# Linker
LD=gcc

PROGRAMS := matrix_test

# all unsere Quelltextdateien und Headerdateien werden in Variablen gespeichert
# dazu wird im momentanen Verzeichnis nach allen .c und .h Dateien gesucht
SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)

# für alle .c-Dateien wird eine Variable erstellt, wo alle .c durch .o ersetzt
# wurden, dies sind unsere Module
MODULES := $(patsubst %.c, %.o, $(SOURCES))


# Objektdateien für Main-module werden separat behandelt
# sie müssen also aus der $(MODULES)-Variable rausgefiltert werden
# dazu nutzen wir filter-out
# die Syntax ist: $(filter-out "zu filternde Objekte", "Objekte aus denen gefiltert werden soll")
# um teilchen.o aus der Liste rauszufiltern, nehmen wir uns $(PROGRAMS) und packen mit
# 'addsuffix' ein .o an jeden Programmnamen
MODULES := $(filter-out $(addsuffix .o, $(PROGRAMS)), $(MODULES))

# benötigte Bibliotheken
LIBS :=

# unser Standardtarget
# wir nutzen es als Gelegenheit, den Inhalt aller Variablen ausgeben zu lassen
all: list-variables $(PROGRAMS)

# das '@' unterdrückt hierbei eine Ausgabezeile, in der die Eingabezeile
# wiederholt wird
list-variables:
	@echo SOURCES=$(SOURCES)
	@echo HEADERS=$(HEADERS)
	@echo MODULES=$(MODULES)
	@echo PROGRAMS=$(PROGRAMS)
	@echo LIBS=$(LIBS)

# Regel zum Kompilieren der Objekdateien in denen die Mainfunktion(en)
# unsere[r,s] Programme sind
# Das/die main-Modul(e) ha[t,ben] keine Headerdateien, haengen aber
# oft von (fast) allen anderen Headerdateien ab
# Die "SECONDEXPANSION" und die doppelten Dollarzeichen sind notwendig, 
# damit die Variablen und Funktionen
# auf der rechten Seite auch wirklich ausgewertet werden koennen...
# https://stackoverflow.com/questions/25611945/getting-the-target-name-in-a-function-in-prerequisites-list-makefile
.SECONDEXPANSION:
$(addsuffix .o, $(PROGRAMS)): % : $$(addsuffix .c,$$(basename $$(@))) $$(HEADERS) Makefile
	$(CC) -c -o $@ $<

# Regel zum kompilieren unserer Module
# Jedes Modul hat eine Headerdatei
%.o: %.c %.h Makefile
	$(CC) -c $<

$(PROGRAMS): % : %.o $(MODULES) Makefile
	$(LD) -o $@ $< $(MODULES) $(LIBS)

.PHONY: clean list-variables
clean:
	rm -f $(MODULES) $(PROGRAMS) $(addsuffix .o, $(PROGRAMS))
