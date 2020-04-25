PROGRAMS := 08_01_quicksort

CC=gcc -ggdb -Wall -Wpedantic -std=c99
LD=gcc
LIBS=-lm

BASENAMES := $(basename $(wildcard *.c))

MODULES := $(filter-out $(PROGRAMS),$(BASENAMES))

all: $(PROGRAMS)

$(PROGRAMS): % :$(addsuffix .o,%) $(addsuffix .o,$(MODULES))
	${LD} -o $@ $^ ${LIBS}

%.o: %.c
	${CC} -c $< -o $@

.PHONY: clean
clean:
	rm $(PROGRAMS) $(addsuffix .o,$(BASENAMES))
