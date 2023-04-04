//------------------------------------------------------------------------------
//-- DUMMY STATEMENT 80 CHARACTERS LONG TO ENSURE THE LINES ARE NOT TOO LONG  --
//------------------------------------------------------------------------------

/**
 * @file graph.c
 * @brief Functions related to undirected graphs
 *
 * @details This module contains functions for manipulating undirected graphs.
 * In particular, there it contains functions for computing the fill-in during
 * a partial or full LU decomposition as well as variants of the minimal degree
 * reordering algorithm.
 *
 * @author Carl Christian Kjelgaard Mikkelsen
 * @version 1.1.0
 * @date 2021-10-6
 * @warning The use of this software is at your own risk
 * @copyright GNU Public Licence
 *
 */

// Libries written primarily for the MD project
#include "graph.h"

#include "list.h"

// Standard libries
#include <stdio.h>
#include <stdlib.h>

void graph_init(graph_t *graph) {
    /**
     * Initializes all pointers in the data structure to NULL
     *
     * @param[in] graph a datastructure representing an undirected graph
     */
    // Nullify pointer to adjacency lists
    graph->adj = NULL;
    // Nullify pointer to row indices
    graph->xadj = NULL;

}

void graph_free(graph_t *graph) {
    /**
     * Releases the heap memory used by all members of the data structure
     *
     * @param[in] graph a datastructure representing an undirected graph
     */
    // Releases the memory used to store a graph.
    free(graph->xadj);
    free(graph->adj);
}

void graph_number_edges(graph_t *graph, int np, int *id,
                        int *perm, int *part) {

    /**
     *  @param[in] graph a pointer to a structure describing a graph G
     *  @param[in] np the number of partitions of the graph
     *  @param[in] id an array such that node k belongs to partition id[k]
     *  @param[out] perm a special permutation, see below
     *  @param[out] part a special partition, see below
     *
     * Given a graph, the edges can be numbered using the given order.
     * This numbering is almost certainly useless for the purpose of doing
     * Gaussian elimination.
     * <br><br>
     * The purpose of this function is to construct a numbering and a partition
     * of the adjacency matrix of the line graph that is suitable for a parallel
     * solve. It is typically possible to refine this numbering using a variant
     * of the minimal degree reordering algorithm.
     * <br><br>
     * It is assumed that id describes a np-way partitioning of G.
     * This partitioning induces a numbering of the edges as follows.
     * The partitions are processed sequentially.
     * First we number all edges that are between vertices in the same partition.
     * Then we number all edges that are between vertices in different partitions.
     * In this way, every edge of G is numbered.
     * <br><br>
     * The array perm is a permutation that describes how to transition from the
     * given ordering of the edges of G to the new numbering.
     * id[k] is the original number of the edge that should be given the number k.
     * part is a an array that describes how to partition the adjacency matrix
     * for the line graph into a arrowhead matrix with k+1 diagonal blocks.
     * The ith (i<k) diagonal block concerns all the edges between vertices in
     * the ith partition of G. The kth diagonal block concerns all edges between
     * vertices of G that are not in the same partition.
     * <br><br>
     * Ideally, the first k diagonal blocks are large and the last block is small.
     * This is what happens when we only need to cut a small number of edges when
     * we partition the graph G. This is not possible in general, but it can be
     * done for a large number of graphs of practical interest.
     *
     * Consider a graph with the following numbered list of edges (a,b)

     * This corresponds to molecule alanine which has m=13 atoms and n=12 bonds.
     * This is a ASCII diagram of the structure of alanine. The atoms/vertices
     * are numbered 0 through 12.
     \verbatim
                7
                |
         0   6--5--8   10-11
          \     |     /
           2----3----9
          /     |     \
         1      4     12
     \endverbatim
     It is very natural to number the bonds as follows
     \verbatim
        id    a    b
       --------------
          0    0    2
          1    1    2
          2    2    3
          3    3    4
          4    3    5
          5    3    9
          6    5    6
          7    5    7
          8    5    8
          9    9   10
         10    9   12
         11   10   11
       --------------
     \endverbatim
     The number is obtained by processing the atoms in increasing order.
     Atom 0 is only bonded to atom 2. Hence bond 0 is between atoms 0 and 2.
     Atom 1 is only bonded to atom 2. Hence bond 1 is between atoms 1 and 2.


     * Now suppose that we do a 2-way partition of the graph/molecule.
     * METIS suggests that the following partitioning
     \verbatim
       Vertex number:   0   1   2   3   4   5   6   7   8   9  10  11  12
       Partition id :   1   1   1   1   1   0   0   0   0   1   1   1   1
     \endverbatim
     * Partition 0 corresponds to the side-chain CH3 and partition 1 is the rest.
     * There is only a single edge between the two partitions, namely a C-C bond.
     * We first number the edges between vertices in partition 0.
     * We then number the edges between vertices in partition 1.
     * Finally, we number the single edge between the two partitions.
     * There is more than one way to achieve this and we simply process the edges
     * in the order they are given above. The sparsity pattern for the adjoint
     * matrix is given below
     \verbatim
         xxx|        |x|
         xxx|        |x|
         xxx|        |x|
         ---------------
            |xxx     | |
            |xxx     | |
            |xxxxx   |x|
            |  xxx   |x|
            |  xxxxx |x|
            |    xxxx| |
            |    xxx | |
            |     x x| |
         ---------------
         xxx|  xxx   |x|
         ---------------
     \endverbatim
     * The 3-by-3 dense block corresponds to the bonds that are internal to the
     * side chain of alanine. The 1-by-1 block in the lower right corner is the
     * single C-C bond that connects the side chain to the carbon back bone.
     * This example is interesting because the fill-in in the second diagonal
     * block is only 2. However, the fill-in in the corresponding off-diagonal
     * blocks is 6. A reordering of the matrix based on a partial minimal degree
     * reordering that respects the partition lines cures the problem. The final
     * matrix has the pattern
     \verbatim
         xxx|        |x|
         xxx|        |x|
         xxx|        |x|
         ---------------
            |x  x    | |
            | xx  x  | |
            | xx  x  | |
            |x  xx  x| |
            |   xx  x| |
            | xx  xxx|x|
            |     xxx|x|
            |   xxxxx|x|
         ---------------
         xxx|     xxx|x|
         ---------------
     \endverbatim
     * While we have lost the pretty overlapping block structure that
     * characterized the initial matrix, the final reordering happens to be
     * perfect and there is no fill-in during the factorization.
     * In terms of the original numbering, the new elimination order is
     \verbatim
     6 7 8 11 0 1 9 10 2 3 5 4
     \endverbatim
     *
     * While alanine is obviously too small for parallel processing, the case
     * of a 4-way partitioning is also interesting. METIS suggests the following
     * partitioning:
     \verbatim
     Vertex number:   0   1   2   3   4   5   6   7   8   9  10  11  12
     Partition id :   3   3   3   1   1   0   0   0   0   1   2   2   1
     \endverbatim
     The sparsity pattern for the adjacency matrix of the line graph is
     \verbatim
         xxx|   | |  | x |
         xxx|   | |  | x |
         xxx|   | |  | x |
         -----------------
            |xx | |  |xx |
            |xxx| |  |xxx|
            | xx| |  |  x|
         -----------------
            |   |x|  |  x|
         -----------------
            |   | |xx|x  |
            |   | |xx|x  |
         -----------------
            |xx | |xx|xx |
         xxx|xx | |  |xx |
            | xx|x|  |  x|
         -----------------
     \endverbatim
     The fill-in during the factorization is 8.
     The second diagonal block does not generate fill-in during a factorization,
     but the two off-diagonal blocks do. The matrix can be reordered as
     \verbatim
     xxx|   | |  | x |
     xxx|   | |  | x |
     xxx|   | |  | x |
     -----------------
        |x x| |  |  x|
        | xx| |  |xx |
        |xxx| |  |xxx|
     -----------------
        |   |x|  |  x|
     -----------------
        |   | |xx|x  |
        |   | |xx|x  |
     -----------------
        | xx| |xx|xx |
     xxx| xx| |  |xx |
        |x x|x|  |  x|
     -----------------
     \endverbatim
     * The fill-in during the factorization is 4.
     * In terms of the original numbering, the new elimination order is
     \verbatim
     6 7 8 10 3 5 11 0 1 2 4 9
     \endverbatim
     *
     * The partitions and sparsity patterns presented in this section where
     * generated used the test program graph_mwe.
     * The relevant call sequences are
     \verbatim
     ./graph_mwe graphs/graph5.txt 2
     ./graph_mwe graphs/graph5.txt 4
     \endverbatim
     *<br><br>
     */

    // For the purpose of the function, we define a proper edge as an edge
    // that is not a self-loop.

    // Isolate the number of vertices
    int m = graph->m;

    // Allocate space for np+1 list
    // Why np+1?
    // The vertices of G have been partitioned into np groups.
    // This induces a partitioning of the line graph L(G)
    // The edges of G fall into two main categories.
    // Vertex a is in group i
    // Vertex b is in group j
    // If i==j, then the edge belongs to partition i of L(G)
    // If i!=j, then the edge belongs to partition np of L(G)

    my_node_t **list = (my_node_t **)malloc((np + 1) * sizeof(my_node_t *));
    // Initialize the head of each list
    for (int i = 0; i < np + 1; i++) {
        list[i] = NULL;
    }

    // Initialize counter for all edges between different vectices
    int s = 0;

    // Loop over the the vertices
    for (int a = 0; a < m; a++) {
        // Loop over the vertices connected to a
        for (int k = graph->xadj[a]; k < graph->xadj[a + 1]; k++) {
            // Isolate the vertex connected to a
            int b = graph->adj[k];
            // Ignore edges from a vertex to itself.
            // WARNING: We assume the graph is undirected.
            if (a < b) {
                // Identify the two relevant partitions
                int i = id[a];
                int j = id[b];
                // Are a and b in the same partition or not?
                if (i == j) {
                    sinsert(s, &list[i]);
                }
                else {
                    sinsert(s, &list[np]);
                }
                // Increment the count of all proper edges
                s++;
            }
        }
    }
    // Initialize an index into perm;
    int k = 0;
    // Now process the lists
    for (int i = 0; i < np + 1; i++) {
        // Compress the ith list into perm starting at position k
        copy_list(list[i], &perm[k]);
        // Update the index pointer
        k += count(list[i]);
    }

    // We can now construct the partitioning of the arrowhead matrix
    part[0] = 0;
    for (int i = 0; i < np + 1; i++) {
        part[i + 1] = part[i] + count(list[i]);
    }

    // Free memory
    for (int i = 0; i < np + 1; i++) {
        free_list(list[i]);
    }
    free(list);
}

//----------------------------------------------------------------------------
void graph_renumber_vertices(graph_t *graph, int *p) {

    /* RENUMBER VERTICES

       Renumbers the vertices of a graph using a given permutation p.
       The permutation sigma is given as in MATLAB

          p = (tau(0), tau(1), tau(2), ...., tau(m-1))

       where tau is the inverse of sigma, i.e. sigma(tau(j)) = j.

       It is necessary to compute the inverse permutation

          q = (sigma(0), sigma(1), sigma(2), ...., sigma(m-1))

       in order to faciliate the transformation.

       REMARK: The MATLAB notation can be annoying; live with it, as there are
       cases where one is preferable to the other and in this case it is simple
       to use both.

    */

    //--------------------------------------------------------------------------
    // Declaration of internal variables
    // -------------------------------------------------------------------------

    // The number of vertices
    int m;

    // Standard counter(s)
    int i, k;

    // Linked list representation of the renumbered graph
    my_node_t **root;

    // The inverse permutation
    int *q;

    // -------------------------------------------------------------------------
    // Start of instructions
    // -------------------------------------------------------------------------

    // Extract the number of vertices in the graph
    m = graph->m;

    // Initialize the linked list representation of the renumbered graph
    root = (my_node_t **)malloc(m * sizeof(graph_t *));
    for (i = 0; i < m; i++) {
        root[i] = NULL;
    }

    // Allocate space for and compute the inverse permutation
    q = (int *)malloc(m * sizeof(int));
    for (i = 0; i < m; i++) {
        q[p[i]] = i;
    }

    // Loop over the adjacency lists and apply the permutation
    for (i = 0; i < m; i++) {
        // Process the ith list
        for (k = graph->xadj[i]; k < graph->xadj[i + 1]; k++) {
            sinsert(q[graph->adj[k]], &root[q[i]]);
        }
    }

    // Release the current graph
    graph_free(graph);

    // Compress the reordered graph into the source
    graph_list_to_compact(m, root, graph);

    free(q);

    // Free the linked lists
    for (i = 0; i < m; i++) {
        free_list(root[i]);
    }
    free(root);

}

// ----------------------------------------------------------------------------
void graph_compute_fill(graph_t *src, graph_t *dest) {
    // Extract the number of vertices from the graph
    int m = src->m;

    // Pointer to the linked list representation of the source graph
    my_node_t **root = (my_node_t **)malloc(m * sizeof(my_node_t *));

    // First expand the graph into a set of linked lists
    graph_compact_to_list(src, root);

    // Loop over the adjacency lists
    for (int i = 0; i < m; i++) {
        // Start at the root of the ith adjancency list
        my_node_t *conductor = root[i];
        // Go over the ith adjacency list one element at a time
        while (conductor != NULL) {
            // We are only interest in higher indices
            if (conductor->number > i) {
                // The apprentice starts at the root of the ith adjancency list
                my_node_t *apprentice = root[i];
                // Let the apprentice run to the end of the list
                while (apprentice != NULL) {
                    // We are only interested in higher indices
                    if (apprentice->number > i) {
                        // Make insertion into list number conductor->number > i
                        sinsert(apprentice->number, &root[conductor->number]);
                        // Please observe that we are NOT changing the ith list !!!
                    }
                    apprentice = apprentice->next;
                }
            }
            conductor = conductor->next;
        }
    }

    // Compress the linked lists into a graph
    graph_list_to_compact(m, root, dest);

    //--------------------------------------------------------------------------
    // Release the dynamically allocated memory
    //--------------------------------------------------------------------------

    // First the m individual lists
    for (int i = 0; i < m; i++) {
        free_list(root[i]);
    }

    // The the root pointer itself
    free(root);
}

void graph_minimal_degree(graph_t *src, int *p) {

    /*

      ALGORITHM

      1: Set Gamma = {0,1,...,m-1} (linked list)
      2: FOR i=0,..,m-1 DO
      3:   Let alpha in Gamma have minimum degree, i.e.

             deg(alpha) = min{deg(gamma)| gamma in Gamma}

      4:   DELETE alpha from alpha's list

      5:   FOR beta IN alpha's list DO
      6:     DELETE alpha from beta's list
      7:   END

      8:   FORM a clique consisting of all nodes in alpha's list

      9:   DELETE alpha from Gamma

     10: END

    */

    //--------------------------------------------------------------------------
    // Declaration of internal variables
    //--------------------------------------------------------------------------

    // Standard counters
    int i;

    // The number of nodes
    int m;

    // The degree of the nodes
    int *deg;

    // Pointers to nodes
    my_node_t **root, *head, *current, *Gamma;

    // A node whose degree is minimal
    int alpha;

    //--------------------------------------------------------------------------
    // Start of instructions
    //--------------------------------------------------------------------------

    // Extract the number of nodes
    m = src->m;

    /* Initialize the list of active nodes

       REMARK: Notice that this need NOT be all available nodes !!!
       Suppose, that we only want a partial factorization and that we want
       nodes a, b, and c to be numbered last. Then we do not add a, b, and c
       to the list and write a, b, and c to the last three entries of the
       permutation.
    */

    /* Build the list of "active" nodes, i.e. nodes which are considered for
       elimination. This list should be shorter if we only want a partial
       factorization */
    Gamma = NULL;
    for (i = 0; i < m; i++) {
        sinsert(i, &Gamma);
    }

    // Allocate space to store the degree of the nodes
    deg = (int *)malloc(m * sizeof(int));

    // Compute the initial degree of the nodes
    for (i = 0; i < m; i++) {
        deg[i] = src->xadj[i + 1] - src->xadj[i];
    }

    //--------------------------------------------------------------------------
    // Build a linked list representation of the graph
    //--------------------------------------------------------------------------

    // Allocate space for the linked list representation of the graph
    root = (my_node_t **)malloc(m * sizeof(my_node_t *));

    // Initialize the individual linked list
    for (i = 0; i < m; i++) {
        root[i] = NULL;
    }

    // Convert the compact representation into linked lists
    graph_compact_to_list(src, root);

    //--------------------------------------------------------------------------
    // Main loop begins here
    //-------------------------------------------------------------------------

    // Initialize the counter.
    i = 0;

    // Loop over the list of active elements until it is empty
    while (Gamma != NULL) {
        // Determine a node alpha which has minimal degree
        alpha = Gamma->number;
        current = Gamma;
        while (current != NULL) {
            if (deg[current->number] < deg[alpha]) {
                alpha = current->number;
            }
            current = current->next;
        }

        // Update the permutation
        p[i] = alpha;

        // Delete alpha from alpha's list
        sdelete(alpha, &root[alpha]);

        // For each beta in alpha's list, delete alpha from beta's list
        current = root[alpha];
        while (current != NULL) {
            deg[current->number] -= sdelete(alpha, &root[current->number]);
            current = current->next;
        }

        // Form a clique of consisting of all members in alpha's list
        head = root[alpha];
        while (head != NULL) {
            current = root[alpha];
            while (current != NULL) {
                deg[current->number] += sinsert(head->number, &root[current->number]);
                current = current->next;
            }
            head = head->next;
        }

        // Delete alpha from the active list
        sdelete(alpha, &Gamma);

        // Increment i in preparation for the next loop.
        i++;
    }

    // Release the adjacency lists
    for (i = 0; i < m; i++) {
        free_list(root[i]);
    }

    // Release the main pointer to the linked list representation of the graph
    free(root);

    // Free the active list
    free_list(Gamma);

    // Release memory for the degrees
    free(deg);
}

//----------------------------------------------------------------------------
void graph_compact_to_list(graph_t *src, my_node_t **root) {

    /* Generates a linked list representation of a graph.

       It is the responsibility of the caller to allocate space for the
       pointers to the head of the lists.

     */

    // Standard counters
    int i, j, k;

    // The number of vertices
    int m = src->m;

    // Initialize the linked lists
    for (i = 0; i < m; i++) {
        root[i] = NULL;
    }

    int count = 0;
    int real_count = 0;
    // Loop over the vertices
    for (i = 0; i < m; i++) {
        // Loop over the elements in the ith adjacency list
        for (k = src->xadj[i]; k < src->xadj[i + 1]; k++) {
            // Which vertex j are we looking at?
            j = src->adj[k];
            // Insert j into the ith list
            real_count += sinsert(j, &root[i]);
            ++count;
        }
    }
}

void graph_list_to_compact(int m, my_node_t **root, graph_t *dest) {

    /* Generates a compact representation of a graph from the linked lists

       DESCRIPTION OF INTERFACE

       ON ENTRY:

       m           the number of vertices/adjacency lists
       root[i]     pointer to the head of the ith list
       dest        pointer to destination datastructure

       ON EXIT
       dest->m     equal to the number of vertices
       dest->xadj
       dest->adj
     */

    // -------------------------------------------------------------------------
    // Declaration of internal variables
    // -------------------------------------------------------------------------

    // Standard counters
    int i;

    // Total length of the adjacency lists AND pointer into dest->adj.
    int length = 0;

    // -------------------------------------------------------------------------
    // Start of instructions
    // -------------------------------------------------------------------------

    // Loop over the lists and compute their total length
    for (i = 0; i < m; i++) {
        length += count(root[i]);
    }

    // Set the number of vertices in the dest
    dest->m = m;

    // Allocate space for dest->xadj;
    dest->xadj = (int *)malloc((m + 1) * sizeof(int));

    // Allocate space for dest->adj
    dest->adj = (int *)malloc(length * sizeof(int));

    // Reset the length/index into dest->xadj
    length = 0;
    dest->xadj[0] = 0;

    // Loop over the lists
    for (i = 0; i < m; i++) {
        // Copy the ith list into dest->adj
        copy_list(root[i], &(dest->adj[length]));
        // Update length
        length += count(root[i]);
        // Fill in the value of dest->xadj[i+1]
        dest->xadj[i + 1] = length;
    }
}
