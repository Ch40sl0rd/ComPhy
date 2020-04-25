SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)
LIBS := -lgsl -lgslcblas

# die erste Regel ist die Standardregel
all: $(basename $(SOURCES))

$(basename $(SOURCES)): % : %.c Makefile
	gcc -g -ggdb -Wall -Wpedantic -std=c99 -o $@ $< ${LIBS}

.PHONY: clean all
clean:
	rm $(basename $(SOURCES))
