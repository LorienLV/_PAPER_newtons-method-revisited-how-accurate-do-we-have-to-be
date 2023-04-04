#include "molecule.h"

#include "list.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

molecule::molecule() {
    invmass = nullptr;
    bonds = nullptr;
    sigmaA = nullptr;
    sigmaB = nullptr;
    sigma2 = nullptr;
    bond_graph = nullptr;
    weights = nullptr;
}

molecule::~molecule() {
    free(invmass);
    free(bonds);
    free(sigmaA);
    free(sigmaB);
    free(sigma2);
    if (bond_graph) {
        graph_free(bond_graph);
        free(bond_graph);
    }
    free(weights);
}

void renumber_bonds(molecule_t *mol, int *p) {

    // -------------------------------------------------------------------------
    // Declaration of internal variables
    // -------------------------------------------------------------------------

    // Number of bonds
    int n;

    // Auxiliary bond list
    int *aux_bonds;
    real *aux_sigmaA;
    real *aux_sigmaB;
    real *aux_sigma2;

    // Standard counter
    int i;

    // -------------------------------------------------------------------------
    // Start of instructions
    // -------------------------------------------------------------------------

    if (mol->bonds != NULL) {
        // Extract the number of bonds
        n = mol->n;

        // Allocate space for the bond auxiliary list
        aux_bonds = (int *)malloc(2 * n * sizeof(int));
        aux_sigmaA = (real *)malloc(n * sizeof(*aux_sigmaA));
        aux_sigmaB = (real *)malloc(n * sizeof(*aux_sigmaB));
        aux_sigma2 = (real *)malloc(n * sizeof(*aux_sigma2));

        // Copy the original bond list into the auxiliary list
        for (i = 0; i < 2 * n; i++) {
            aux_bonds[i] = mol->bonds[i];
        }
        for (i = 0; i < n; i++) {
            aux_sigmaA[i] = mol->sigmaA[i];
            aux_sigmaB[i] = mol->sigmaB[i];
            aux_sigma2[i] = mol->sigma2[i];
        }

        // Apply the permutation to the bond list
        for (i = 0; i < n; i++) {
            mol->bonds[2 * i + 0] = aux_bonds[2 * p[i] + 0];
            mol->bonds[2 * i + 1] = aux_bonds[2 * p[i] + 1];

            mol->sigmaA[i] = aux_sigmaA[p[i]];
            mol->sigmaB[i] = aux_sigmaB[p[i]];
            mol->sigma2[i] = aux_sigma2[p[i]];
        }

        // Free the auxiliary bond list
        free(aux_bonds);
        free(aux_sigmaA);
        free(aux_sigmaB);
        free(aux_sigma2);
    }
    if (mol->bond_graph != NULL) {
        // Apply the permutation to the bond graph
        graph_renumber_vertices(mol->bond_graph, p);
    }
}

//----------------------------------------------------------------------------
void make_bond_graph(molecule_t *mol) {
    /* Builds the bond graph from the bond list

       Description of interface

       ON ENTRY:
           mol->m       the number of atoms
           mol->n       the number of bonds
           mol->bonds   sequential list of n bonds

       ON EXIT:
           mol->graph   a pointer to a compact representation of the bond graph

       No other variables are modified!

       In the bond graph there is an edge between vertices i and j if and only
       if bonds i and j have one atom in common. The bond graph is stored using
       two arrays called XADJ and ADJ. Each adjacency list is stored in strictly
       increasing order. The n lists are stored in strictly increasing order in
       the array ADJ. The number XADJ[j] is the index inside ADJ of the first
       entry of the jth adjacency list. XADJ[n] points just beyond the end of
       ADJ and is the length of ADJ.

       IDEA:

       Every atom participates in at most K bonds. Suppose that we have m
       auxiliary lists, such that the ath list gives the bonds which atom a
       partakes in. Then for each bond i, involving atoms a(i) and b(i), we can
       immediately construct the adjacency list for bond i, simply by merging
       all the auxiliary list for atoms a(i) and b(i). Therefore our first step
       is the construction of these auxiliary lists.

       ALGORITHM

       STEP 1: Construction of the auxiliary lists

       for each bond i do
         for each atom a in bond i do
           insert i into the ath auxiliary list
         end
       end

       COST: There are n bonds which involve 2 atoms each. Each auxiliary list
       will never be longer than K, and inserting an element into a sorted list
       of length L requires at most L comparisons, so the cost is less than 2nK
       comparisons.

       STEP 2: Construction of the adjacency lists

       for each bond i do
         for each atom a in bond i do
           for each bond j in the ath auxiliary list do
             insert j into the ith adjacency lists
           end
         end
       end

       COST: There are again n bonds which involve at 2 atoms each. Each
       auxiliary list has a length which is at most K. Each adjacency list will
       have a length which is at most 2K+1. Inserting M elements into a list of
       length at most L requires less than LM comparisons, so in total we do
       O(nK^2) comparisons.

       STEPS 3: Compress the adjacency lists into the array ADJ and create XADJ.
       There are n adjacency lists of length at most 2K+1, so the cost is O(nK).

    */


    // -------------------------------------------------------------------------
    // Declaration of internal variables
    // -------------------------------------------------------------------------

    // The number of atoms
    int m;

    // The number of bonds
    int n;

    // The list of bonds
    int *bonds;

    // The bond graph
    graph_t *graph;

    // Standard counters
    int i, j;

    // The current bond is between atoms a and b.
    int a, b;

    // An auxilary variable used to generate XADJ
    int temp;

    // Variables needed to manipulate linked lists
    struct my_node *conductor;
    struct my_node **aux, **root;

    // -------------------------------------------------------------------------
    // Start of instructions
    // -------------------------------------------------------------------------

    // Extract the number of atoms (m) and the number of bonds
    m = mol->m;
    n = mol->n;

    // Establish shortcut to the list of bonds
    bonds = mol->bonds;

    // Allocate space for the bond graph
    graph = (graph_t *)malloc(sizeof(graph_t));

    // allocate space for the auxilary lists
    aux = (struct my_node **)malloc(m * sizeof(struct my_node *));

    /* Initialize the auxiliary lists. The ith auxiliary list will record the
       bonds which atom i partakes in. */

    for (i = 0; i < m; i++) {
        aux[i] = NULL;
    }

    // Loop over the sequential list of bonds
    for (j = 0; j < n; j++) {
        // Isolate the numbers of the atoms which partake in bond j.
        a = bonds[2 * j];
        b = bonds[2 * j + 1];
        // Insert bond j into the ath auxiliary list
        sinsert(j, &aux[a]);
        // Insert bond j into the bth auxiliary list
        sinsert(j, &aux[b]);
    }

    /* This completes the construction of the auxiliary list. For each atom we
       now have a list of the bonds it partakes in! */

    // Allocate space for n adjacency lists
    root = (struct my_node **)malloc(n * sizeof(struct my_node *));

    // Initialize the adjacency lists
    for (j = 0; j < n; j++) {
        root[j] = NULL;
    }

    // Loop over the sequential list of bonds
    for (j = 0; j < n; j++) {
        // Isolate the numbers of the atoms which partake in bond j.
        a = bonds[2 * j];
        b = bonds[2 * j + 1];
        /* Insert every element of the a'th auxiliary list into the adjacency
           list for bond j. */
        conductor = aux[a];
        while (conductor != NULL) {
            sinsert(conductor->number, &root[j]);
            conductor = conductor->next;
        }
        /* Insert every element of the b'th auxiliary list into the adjacency
           list for bond j */
        conductor = aux[b];
        while (conductor != NULL) {
            sinsert(conductor->number, &root[j]);
            conductor = conductor->next;
        }
    }

    // Free the memory used by the auxilliary lists
    for (i = 0; i < m; i++) {
        free_list(aux[i]);
    }
    // Do not forget to release aux itself.
    free(aux);

    /* This complete the construction of the adjacency lists as simply linked
       list. It remains to compress the information into two arrays XADJ and
       ADJ. */

    // Allocate space for the array of indices into the combined adjacency list
    graph->xadj = (int *)malloc((n + 1) * sizeof(int));

    // Initialize the counters
    temp = 0;
    graph->xadj[0] = temp;

    // Find the length of the adjacency lists and compute the entries of XADJ
    for (j = 0; j < n; j++) {
        temp = temp + count(root[j]);
        graph->xadj[j + 1] = temp;
    }

    // Allocate space for the graph
    graph->adj = (int *)malloc(temp * sizeof(int));

    // Copy the lists into graph->adj and free the memory used.
    for (j = 0; j < n; j++) {
        copy_list(root[j], &graph->adj[graph->xadj[j]]);
        free_list(root[j]);
    }
    // Do not forget to release root itself
    free(root);

    // Set the number of vertices in the bond graph
    graph->m = n;

    // Do NOT forget to save the result of your work into MOL!
    mol->bond_graph = graph;
}

//----------------------------------------------------------------------------
void make_weights(molecule_t *mol, graph_t *graph) {

    /* Precomputes the weights needed to generate the matrix A(x,y).

       Description of interface:

       ON ENTRY:

         mol->m        the number of atoms
         mol->invmass  pointer to the inverse of the atomic masses
         mol->n        the number of bonds
         mol->bonds    a pointer to a sequential list of n bonds
         graph         a pointer to a compact representation of the graph
                       to use, typically the bond graph or the lower triangular
                       portion of it

      ON EXIT:
         mol->weights  a pointer to list of weights compatible with the graph

      No other variables have been modified.

      REMARK(S):

       1) At this point we do not take fill into account. Fill is identified
          by playing the elimination game on the bond graph. Right now, I am
          simply not sure where this code should go.


       Description of the construction of the weights:

       The following example explains the origins of the weights and how to
       compute them.

       EXAMPLE:

       A molecule with m = 6 atoms and n = 5 bonds. The atoms are numbered 0
       through 5, and the bonds are numbered 0 through 4. It is irrelevant if
       this molecule exists in the real world or not :)


        0                 4           The bonds are given by the array
         \               /
         (0)           (3)            bonds = {0,2; 1,2; 2,3; 3,4; 3,5}
           \           /
            2 ------- 3               Each bond involves 2 atoms. Bond 2 is
           /    (2)    \              between atoms 2 and 3.
         (1)           (4)
         /               \
       1                  5

       The adjacency graph for the bonds is the graph


          0   3              matrix A is 5 by 5
         / \ / \
        |   2   |            |xxx  |
         \ / \ /             |xxx  |
          1   4              |xxxxx|
                             |  xxx|
                             |  xxx|

       This is how the graph is encoded

        adj  = {0,1,2; 0,1,2; 0,1,2,3,4; 2,3,4; 2,3,4}     (17 entries)
       xadj  = {0, 3, 6, 11, 14, 17}                       ( 6 entries)

       Please note that there are n+1 entries in XADJ, and the last entry point
       just beyond the end of adj. Therefore XADJ[n+1] is the number of nonzero
       entries in the matrix.

       NOTATION:

         1) We write rij or r(i,j) for the vector from atom i to j.
         2) We write <x,y> for the scalar product between the vectors x and y.

       In the above example bond number 2 involves atoms 2 and 3 and is a
       (mathematical) constraint of the type

                  0.5*(||r23||^2 - (bond length)^2) = 0

       The factor 0.5 is only included to give the Jacobian of the constraint
       function a form which I (CCKM) consider aesthetically pleasing.

       Now, our matrix of the form

                  A =  Dg(r)*inv(M)*Dg(s)'

       where r and s are vectors describing two different configurations of
       the m atoms, so r and s each have 3m components each.

       Below is a table of the nonzero entries of the matrix A for our current
       example:

                   bond = {0,2,1,2,2,3,3,4,3,5}

       i    j      entry A(i,j)
       ---------------------------------------------------------------
       0    0     (invmass(0)+invmass(2))*<r02,s03>
       1    0     +invmass(2)*<r12,s02>
       2    0     -invmass(2)*<r23,s02>

       0    1     +invmass(2)*<r02,s12>
       1    1     (invmass(1)+invmass(2))*<r12,s12>
       2    1     -invmass(2)*<r23,s12>

       0    2     -invmass(2)*<r02,s23>
       1    2     -invmass(2)*<r12,s23>
       2    2     (invmass(2)+invmass(3))*<r23,s23>
       3    2     -invmass(3)*<r34,s23>
       4    2     -invmass(3)*<r35,s23>

       2    3     -invmass(3)*<r23,s34>
       3    3     (invmass(3)+invmass(4))*<r34,s34>
       4    3     +invmass(3)*<r35,s34>

       2    4     -invmass(3)*<r23,s35>
       3    4     +invmass(3)*<r34,s35>
       4    4     (+invmass(3)+invmass(5))*<r35,s35>

       This table is very carefully constructed! Please note the following:

       a) Reading the table from the top to the bottom takes us through the
       nonzero entries of A in column major format, exactly as matrices are
       stored in LAPACK. Moreover, we are writing to main memory in an order
       which is strictly increasing.

       b) For each column of the matrix A we deal with a FIXED vector s(a,b),
       rather than both s(a,b) and s(b,a) = -s(a,b).

       c) The order of the indices as in r(a,b) for a bond k, is EXACTLY the
       order in which the atoms which partake in bond k are given in the bond
       table.

       EXAMPLE: As noted above bond 2 is a bond between atoms 2 and 3. The bond
       table lists atom 2 BEFORE atom 3. Therefore, we use the vector r23,
       rather than the (equvalent) vector r32 = -r23

       d) The diagonal entries are different from the off diagonal entries
       because the weights are different.

       In order to efficiently generate the matrix A we precompute the following
       weights

       i    j     weight
       ----------------------------------------------
       0    0     +invmass(0)+invmass(2)
       1    0     +invmass(2)
       2    0     -invmass(2)

       0    1     +invmass(2)
       1    1     +invmass(1)+invmass(2)
       2    1     -invmass(2)

       0    2     -invmass(2)
       1    2     -invmass(2)
       2    2     +invmass(2)+invmass(3)
       3    2     -invmass(3)
       4    2     -invmass(3)

       2    3     -invmass(3)
       3    3     +invmass(3)+invmass(4)
       4    3     +invmass(3)

       2    4     -invmass(3)
       3    4     +invmass(3)
       4    4     +invmass(3)+invmass(5)

       In short, the weights includes all the relevant information about the
       signs as well as the masses of the atoms.

       Given vectors r and s with a 3m components as well as the bond graph, we
       can now generate the matrix A one entry at a time, moving through RAM
       memory in a strictly increasing order.

    */

    // -------------------------------------------------------------------------
    // Declaration of internal variables
    // -------------------------------------------------------------------------

    // The number of bonds;
    int n;

    // Shortcut to the list of inverse masses
    real *invmass;

    // Shortcut to the list of bonds
    int *bonds;

    // Shortcut to the list of weights
    real *weights;

    // The number of nonzeros
    int nnz;

    // Standard counters
    int i, j, k;

    // Atomic labels
    int a, b, c, d, e;

    // -------------------------------------------------------------------------
    // Start of instructions
    // -------------------------------------------------------------------------

    // Extract the number of bonds
    n = mol->n;

    // Establish shortcut to the list of inverse masses.
    invmass = mol->invmass;

    // Establish shortcut to the list of bonds.
    bonds = mol->bonds;

    // Extract the number of nonzero entries
    nnz = graph->xadj[n];

    // Allocate enough space for nnz weights
    weights = (real *)malloc(nnz * sizeof(weights));

    // Initialize the pointer to the nonzero entries
    k = 0;
    // Loop over the bonds or equivalently the ROWS of the matrix
    for (i = 0; i < n; i++) {
        // Isolate the atoms which partake in the ith bond
        a = bonds[2 * i];
        b = bonds[2 * i + 1];
        // Loop over the entries of the adjacency list of the ith bond.
        for (j = graph->xadj[i]; j < graph->xadj[i + 1]; j++) {
            if (i != graph->adj[j]) {
                // This is an off diagonal diagonal entry.
                // We begin by isolating the atoms which partake in the ith bond
                c = bonds[2 * graph->adj[j]];
                d = bonds[2 * graph->adj[j] + 1];

                // Is a fillin.
                if (a != c && a != d && b != c && b != d) {
                    weights[k] = 0;
                }
                else {
                    // Determine the atom which is common to bond j and bond i
                    if ((a == c) || (a == d)) {
                        // The common atom is atom a
                        e = a;
                    }
                    else {
                        // The common atom must necessarily be atom b.
                        e = b;
                    }
                    weights[k] = invmass[e];
                    // Determine the appropriate sign
                    if ((a == d) || (b == c)) {
                        /* You should reverse the order the atoms for one of the two bonds,
                        but this impractical, so we just reverse the sign of the weight
                        */
                        weights[k] = -weights[k];
                    }
                }
            }
            else {
                /* This is a diagonal entry.
                   The weight is special, but the order of the atoms in the bond list
                   is irrelevant. Yes, you could flip the order of one pair of atoms,
                   but then you would be compelled to flip the order of the second
                   pair, and so you would change sign twice.
                */
                weights[k] = invmass[a] + invmass[b];
            }

            // Move on to the next nonzero entry of the matrix
            k++;
        }
    }
    // Remember to save the results of your work!
    mol->weights = weights;
}
