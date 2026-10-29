# SIS simulation example

Adapted from Rui Costa's `SIS_localization_l_0.12_N_10^6_runs_50000.c` (2016).
The program simulates an SIS process. Extend it to report an infection
occupation probability for every vertex over the retained runs.
The program implements the standard SIS model with recovery parameter 1
and infection parameter lambda per neighbour.

Estimate the fraction of process time each vertex is infected under the standard SIS dynamics.

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
of initial runs excluded from observation (0..runs-1), and seed (1..2147483647).
The program returns nonzero if its iteration limit is reached.

Replay arguments are the graph, lambda and trace path.
Records contain the sequence of observed configurations.

## Output contract

In simulation mode, the existing stderr contains a single summary line with vertices, runs,
burnin and iterations. The `count_hist COUNT VALUE` rows report the fraction
of iterations in retained runs that begin with COUNT infected vertices,
for COUNT=1..N.

Add one stdout row per vertex:

`vertex ID PROBABILITY`

Vertex labels remain fixed; row order is unspecified. Preserve the existing
commands, seeded simulation behavior, trace output, count_hist rows and summary. The
extension must not introduce additional random draws. Both simulation and
replay must report the new observable. Successful execution returns zero;
invalid input or an incomplete run returns nonzero. Supply source and build
instructions.
