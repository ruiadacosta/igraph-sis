#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <igraph.h>

#define MAX_N 128
#define POOL_SIZE 64
#define MAX_ITERATIONS 5000000ULL

typedef struct {
    double lambda;
} ModelConfig;

static ModelConfig model;

typedef struct {
    int n, runs, burnin;
    unsigned long long count_hist[MAX_N + 1];
    unsigned long long norm;
    FILE *trace;
} Stats;

static void fail(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(2);
}
static uint64_t rng_state;
static double uniform(void) {
    uint64_t x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return (double)((x * UINT64_C(2685821657736338717)) >> 11) *
           (1.0 / 9007199254740992.0);
}
static int choose(int n) { return (int)(uniform() * n); }

static void record_state(Stats *s, int run, const int *active, int count) {
    if (!count)
        return;

    if (s->trace) {
        fprintf(s->trace, "%d %d", run, count);
        for (int j = 0; j < count; ++j) fprintf(s->trace, " %d", active[j]);
        fputc('\n', s->trace);
    }
    if (run >= s->burnin) {
        ++s->count_hist[count];
        ++s->norm;
    }
}
static void report(const Stats *s) {
    if (!s->norm) fail("No records to report");
    for (int k = 1; k <= s->n; ++k)
        printf("count_hist %d %.17g\n", k, (double)s->count_hist[k] / (double)s->norm);
    fprintf(stderr, "vertices=%d runs=%d burnin=%d records=%llu\n",
            s->n, s->runs, s->burnin, s->norm);
}

static long parse_integer(const char *text, long lo, long hi) {
    char *end;
    errno = 0;
    long v = strtol(text, &end, 10);
    if (errno || !*text || *end || v < lo || v > hi) fail("Invalid integer argument");
    return v;
}
static double parse_lambda(const char *text) {
    char *end;
    errno = 0;
    double v = strtod(text, &end);
    if (errno || !*text || *end || !isfinite(v) || v < 0.0 || v > 10.0)
        fail("Invalid lambda");
    return v;
}
static long next_integer(char **cursor) {
    char *p = *cursor, *end;
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;
    errno = 0;
    long v = strtol(p, &end, 10);
    if (errno || end == p) fail("Invalid input row");
    *cursor = end;
    return v;
}
static void end_row(char *p) {
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;
    if (*p) fail("Extra input columns");
}
typedef struct { int source, target; } Attempt;
static int n;
static igraph_t graph;
static igraph_adjlist_t adjacency;
static int active[MAX_N], infected[MAX_N];
static int infected_count, attempt_count;
static Attempt attempts[MAX_N * (MAX_N - 1)];

static int vertex_degree(int v) {
    return (int)igraph_vector_int_size(igraph_adjlist_get(&adjacency, v));
}
static int vertex_neighbor(int v, int j) {
    return (int)VECTOR(*igraph_adjlist_get(&adjacency, v))[j];
}

static void read_graph(const char *path) {
    igraph_vector_int_t edges;
    if (igraph_vector_int_init(&edges, 0)) fail("Require a simple undirected graph");
    n = 0;
    FILE *f = fopen(path, "r");
    if (!f) fail("Cannot open graph");
    char line[4096];
    while (fgets(line, sizeof line, f)) {
        char *p = line;
        long u = next_integer(&p), v = next_integer(&p);
        end_row(p);
        if (u < 0 || u >= MAX_N || v < 0 || v >= MAX_N || u == v)
            fail("Require a simple undirected graph");
        if (igraph_vector_int_push_back(&edges, u) ||
            igraph_vector_int_push_back(&edges, v)) fail("Require a simple undirected graph");
        if (u >= n) n = (int)u + 1;
        if (v >= n) n = (int)v + 1;
    }
    if (ferror(f)) fail("Graph read error");
    fclose(f);
    if (n < 3 || igraph_vector_int_size(&edges) == 0)
        fail("Require 3..128 vertices");
    if (igraph_create(&graph, &edges, n, IGRAPH_UNDIRECTED))
        fail("Require a simple undirected graph");
    igraph_vector_int_destroy(&edges);
    igraph_bool_t simple, connected;
    if (igraph_is_simple(&graph, &simple, IGRAPH_UNDIRECTED) || !simple)
        fail("Require a simple undirected graph");
    if (igraph_is_connected(&graph, &connected, IGRAPH_WEAK) || !connected)
        fail("Require connected graph and dense vertex labels");
    if (igraph_adjlist_init(&graph, &adjacency, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE))
        fail("Require a simple undirected graph");
}

static void add_active_vertex(int v) {
    active[infected_count++] = v;
    for (int j = 0; j < vertex_degree(v); ++j)
        attempts[attempt_count++] = (Attempt){v, vertex_neighbor(v, j)};
}

static void remove_active_vertex(int v) {
    int j = 0;
    int begin = 0;
    while (active[j] != v) {
        begin += vertex_degree(active[j]);
        ++j;
    }

    int last = active[infected_count - 1];
    if (j != infected_count - 1) {
        int tail = attempt_count - vertex_degree(last);
        int middle = tail - begin - vertex_degree(v);
        Attempt saved[MAX_N];

        memcpy(saved, attempts + tail, (size_t)vertex_degree(last) * sizeof *saved);
        memmove(attempts + begin + vertex_degree(last),
                attempts + begin + vertex_degree(v),
                (size_t)middle * sizeof *attempts);
        memcpy(attempts + begin, saved, (size_t)vertex_degree(last) * sizeof *saved);
    }

    attempt_count -= vertex_degree(v);
    active[j] = active[--infected_count];
}

static int recovery_selected(void) {
    double ratio =
        model.lambda * ((double)attempt_count / infected_count);

    return uniform() < 1.0 / (1.0 + ratio);
}

static void process_attempt(Stats *s, int run) {
    int v, value;
    if (recovery_selected()) {
        v = active[choose(infected_count)];
        value = 0;
    } else {
        v = attempts[choose(attempt_count)].target;
        value = 1;
    }

    value ^= infected[v];
    infected[v] ^= value;

    record_state(s, run, active, infected_count * value);

    if (value) {
        if (infected[v])
            add_active_vertex(v);
        else
            remove_active_vertex(v);
    }
}

static void replay(const char *graph_path, const char *lambda_text, const char *path) {
    read_graph(graph_path);
    model.lambda = parse_lambda(lambda_text);
    FILE *f = fopen(path, "r");
    if (!f) fail("Cannot open trace");
    Stats s = {0};
    char line[4096], extra;
    if (!fgets(line, sizeof line, f) ||
        sscanf(line, "sis_trace vertices %d runs %d burnin %d %c", &s.n, &s.runs,
               &s.burnin, &extra) != 3 || s.n != n ||
        s.runs < 1 || s.runs > 100000 || s.burnin < 0 || s.burnin >= s.runs)
        fail("Invalid trace header");
    int previous = -1;
    while (fgets(line, sizeof line, f)) {
        if (!strchr(line, '\n') && !feof(f)) fail("Trace row too long");
        char *p = line;
        long run = next_integer(&p);
        long count = next_integer(&p);
        if (run < 0 || run >= s.runs || run < previous || run > previous + 1 ||
            count < 1 || count > s.n) fail("Invalid trace state");
        int active[MAX_N], seen[MAX_N] = {0};
        for (int j = 0; j < count; ++j) {
            long v = next_integer(&p);
            if (v < 0 || v >= s.n || seen[v]) fail("Invalid active vertex");
            active[j] = (int)v;
            seen[v] = 1;
        }
        end_row(p);
        previous = (int)run;
        record_state(&s, (int)run, active, (int)count);
    }
    if (ferror(f) || previous != s.runs - 1) fail("Incomplete trace");
    fclose(f);
    if (!s.norm) fail("No records to report");
}

static void simulate(int argc, char **argv) {
    if (argc != 6 && argc != 8) fail("Usage: sis_sim graph lambda runs burnin seed [--trace path]");
    int runs = (int)parse_integer(argv[3], 1, 100000);
    int burnin = (int)parse_integer(argv[4], 0, runs - 1);
    rng_state = (uint64_t)parse_integer(argv[5], 1, 2147483647L);
    model.lambda = parse_lambda(argv[2]);
    read_graph(argv[1]);
    Stats s = {0};
    s.n = n; s.runs = runs; s.burnin = burnin;
    if (argc == 8) {
        if (strcmp(argv[6], "--trace")) fail("Unknown option");
        if (!strcmp(argv[7], argv[1])) fail("Trace path must differ from graph path");
        s.trace = fopen(argv[7], "w");
        if (!s.trace) fail("Cannot create trace");
        fprintf(s.trace, "sis_trace vertices %d runs %d burnin %d\n", n, runs, burnin);
    }
    int pool[POOL_SIZE][MAX_N], lengths[POOL_SIZE];
    for (int i = 0; i < POOL_SIZE; ++i) { pool[i][0] = n - 1; lengths[i] = 1; }
    unsigned long long i = 0;
    for (int run = 0; run < runs; ++run) {
        int slot = choose(POOL_SIZE);
        infected_count = attempt_count = 0;
        memset(infected, 0, sizeof infected);
        for (int j = 0; j < lengths[slot]; ++j) {
            int v = pool[slot][j];
            infected[v] = 1;
            add_active_vertex(v);
        }
        while (infected_count) {
            if (++i > MAX_ITERATIONS) fail("Simulation event limit reached");
            if (uniform() < 1.0 / POOL_SIZE) {
                slot = choose(POOL_SIZE);
                lengths[slot] = infected_count;
                memcpy(pool[slot], active, (size_t)infected_count * sizeof(int));
            }
            process_attempt(&s, run);
        }
    }
    if (s.trace && fclose(s.trace)) fail("Trace write error");
    report(&s);
}
int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "--replay")) {
        if (argc != 5) fail("Usage: sis_sim --replay graph lambda trace");
        replay(argv[2], argv[3], argv[4]);
    } else simulate(argc, argv);
    igraph_adjlist_destroy(&adjacency);
    igraph_destroy(&graph);
    return 0;
}
