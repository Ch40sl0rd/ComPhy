SOURCES := $(wildcard *.c)

all: $(basename $(SOURCES))

$(basename $(SOURCES)): % : %.c Makefile
	gcc -g -ggdb -Wall -Wpedantic -std=c99 -o $@ $< -lm

.PHONY: clean
clean:
	rm $(basename $(SOURCES))
