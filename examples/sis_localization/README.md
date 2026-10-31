# SIS simulation example

The program simulates an SIS process. Estimate each vertex’s infection
probability under the standard SIS dynamics.
The program implements the standard SIS model with recovery parameter 1
and infection parameter lambda per neighbour.

## Build and run

Requires a C99 compiler, make, pkg-config and the igraph C library.

```
make
./sis_sim fixtures/sample_a.edgelist 0.12 300 30 111
./sis_sim fixtures/sample_a.edgelist 0.12 300 30 111 --trace trajectory.tsv
./sis_sim --replay fixtures/sample_a.edgelist 0.12 trajectory.tsv
```

Graph inputs contain one undirected edge per line, with dense vertex labels
0 through N-1. They must be simple and connected with 3..128 vertices.
Simulation arguments are lambda (0..10), number of runs (1..100000), number
of initial runs excluded from the reported results (0..runs-1), and seed (1..2147483647).
The program returns nonzero if its iteration limit is reached.

Replay arguments are the graph, lambda and trace path.
The trace contains a sequence of recorded configurations.

## Output contract

In simulation mode, stderr contains a single summary line with vertices, runs,
burnin and records. The `count_hist COUNT VALUE` rows report the fraction of
reported records containing COUNT infected vertices, for COUNT=1..N.

Replay reports the requested vertex profiles only; it does not print
`count_hist` rows or the simulation summary.

Add one stdout row per vertex:

`vertex ID PROBABILITY`

Vertex labels remain fixed; row order is unspecified. Preserve the existing
commands, trace format, seeded trajectories and random-number call sequence.
Follow the histogram and summary definitions above. Both simulation and replay
must report the new observable. Successful execution returns zero; invalid
input or an incomplete run returns nonzero. Supply source and build
instructions.
