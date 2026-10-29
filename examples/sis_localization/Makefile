CC ?= cc
CFLAGS ?= -std=c99 -O2 -Wall -Wextra -Werror
IGRAPH_CFLAGS := $(shell pkg-config --cflags igraph)
IGRAPH_LIBS := $(shell pkg-config --libs igraph)

all: sis_sim
sis_sim: sis_sim.c
	$(CC) $(CFLAGS) $(IGRAPH_CFLAGS) sis_sim.c -o sis_sim $(IGRAPH_LIBS) -lm
clean:
	rm -f sis_sim
.PHONY: all clean
