//----------------------------------------------------------------------------
//-- DUMMY STATEMENT 78 CHARACTERS LONG TO ENSURE THE LINES ARE NOT TOO LONG -
//----------------------------------------------------------------------------

/*  GRAPH, version 1.0.0

    Author:
    Carl Christian Kjelgaard Mikkelsen,
    Department of Computing Science and HPC2N
    Umeaa University
    Sweden
    Email: spock@cs.umu.se

    Date: August 14th, 2014

    This is a module for manipulating general graphs (vertices and edges).

*/
#ifndef GRAPH_GUARD
#define GRAPH_GUARD

// TODO: THIS SHOULD BE C++. Graph should be a class or atleast a C++ struct.

// Project specific modules
#include "list.h"

// ---------------------------------------------------------------------------
// Definition of datastructures
// ---------------------------------------------------------------------------

// Data structure for graphs
typedef struct graph {
    /*
     * Here m is the number of vertices, and xadj[i] is the index in adj of the
     * first element in the ith adjency list. The adjacency lists are stored
     * contiguously in the array adj. This is a very compact way of storing
     * graphs.
     *
     * Example:
     *
     * 0 - 1
     * |
     * 2
     *
     * m = 3
     * xadj = [0, 3, 5, 7]
     * adj = [0, 1, 2, 0, 1, 0, 2]
     *
     */
    int m;
    int *xadj;
    int *adj;
} graph_t;

/* HISTORICAL REMARK:

   Previous I was not familiar with the use of the TYPEDEF keyword, so a few
   words are in order. This constructions allows us to issue a declaration of
   the form

      graph_t   my_graph;

   and we can now refer to the variables my_graph.m, my_graph.xadj, and
   my_graph.adj. Without the use of the TYPEDEF keyword, we would have to
   write

      struct graph  my_graph;

   in order to acheive the same effect.

*/

// ---------------------------------------------------------------------------
// Function prototypes
// ---------------------------------------------------------------------------

// TENTATIVE: Build a large graph from a set of subgraphs
// void build_graph(graph_t **graph, graph_t **connections);

// Initialize a graph
void graph_init(graph_t *graph);

// Releases the memory held by a given graph
void graph_free(graph_t *graph);

// Minimal degree ordering where equivalent nodes are eliminated together.
void graph_minimal_degree(graph_t *src, int *p);

// Renumber the nodes of a graph using a given permutation
void graph_renumber_vertices(graph_t *graph, int *p);

// Computes the fill-in graph from the source graph
void graph_compute_fill(graph_t *src, graph_t *dest);

// Makes a set of linked lists from a compressed representation of a graph
void graph_compact_to_list(graph_t *src, my_node_t **root);

// Makes a compressed representation of a graph from a set of linked lists
void graph_list_to_compact(int m, my_node_t **root, graph_t *dest);

#endif
