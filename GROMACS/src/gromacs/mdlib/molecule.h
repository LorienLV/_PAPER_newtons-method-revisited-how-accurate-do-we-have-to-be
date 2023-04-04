#pragma once

#include "gromacs/utility/real.h"

#include "graph.h"

/* MOLECULE

	 Molecules structures and routines.

	 Main author:
	 Carl Christian Kjelgaard Mikkelsen,
	 Department of Computing Science and HPC2N
	 Umeaa University
	 Sweden
	 Email: spock@cs.umu.se

	 Additional programming by:
	 Jesús Alastruey Benedé
	 Departamento de Informática e Ingeniería de Sistemas
	 Universidad de Zaragoza
	 España
	 Email: jalastru@unizar.es

	 Additional programming by:
	 Lorién López Villellas
	 Departamento de Informática e Ingeniería de Sistemas
	 Universidad de Zaragoza
	 España
	 Email: lorien.lopez@unizar.es

*/

/* This structure contains all the information which can be derived from the
	 bond list and the list of masses. As soon as x, and y are given, this is
	 the information you need to compute the matrix A = A(x,y) rapidly with
	 monotone access of the main memory.

	 You will only need one structure for each type of molecule present in your
	 simulation.
*/
typedef struct molecule {
    int m;                  // The number of atoms
    int n;                  // number of bonds
    real *invmass;          // invmass[i] is 1/mass of the ith atom
    int *bonds;             // sequential bond list similar to GROMACS

    real *sigmaA;            // sequential list of bond lengths A
    real *sigmaB;            // sequential list of bond lengths B
    real *sigma2;            // sequential list of the square of bond lengths A
    //
    graph_t *bond_graph;    // the line graph of the atomic graph
    real *weights;          // weights used to compute the entries of A(x,y)

    /**
     * A constructor that sets every pointer to null, to make this struct
     * easier to manage.
     */
    molecule();

    /**
     * Free dynamically allocated data.
     *
     */
    ~molecule();
} molecule_t;

//----------------------------------------------------------------------------
/* Renumber the bonds (list/graph) using a given permutation.

    There is a important point to consider here:

    Do we or do we not allow subroutines to change pointer values when we
    are simply updating the information which they target?

    Here it would have been easy to save a 2*n memory operations by simply
    NOT copying the bond list into the auxiliary list and simply changing
    the value of the pointer mol->bond_graph!

    However, CCKM finds that it is far safer not to change this pointer.
    Suppose the user had created another pointer to the bond graph before
    calling renumber_bonds. This pointer **might** be invalid on the return
    from renumber_bonds.

    Tentative guiding principles:

        1) If the target size is unchanged, then preserve the pointer
        2) If the target size might change, then do not preserve the pointer.

    WARNING: There is only limited check of the sanity of the input.
*/
void renumber_bonds(molecule_t *mol, int *p);

// Builds the bond graph from the bond list
void make_bond_graph(molecule_t *mol);

// Computes auxiliary weights needed for the construction of A
void make_weights(molecule_t *mol, graph_t *graph);
