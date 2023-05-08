#pragma once

#include <fstream>
#include <vector>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/real.h"
#include "molecule.h"

/**
 * @author Lorién López-Villellas (lorien.lopez@unizar.es)
 *
 */

namespace gmx {

class Ilves {
public:
    /**
     * Initializes the ILVES solver.
     *
     * @param mol Molecule structure.
     */
    Ilves(molecule_t *mol);

    /**
     * Destructor.
     *
     */
    ~Ilves();

    /**
     * Try to solve the bond constraints of MOL in at most MAXIT iterations of
     * Newton's method with a tolerance for each atom of TOL.
     *
     * @param time_step The time step.
     * @param niters When returning, NITERS will contain the total number of
     * iterations executed by the solver.
     *
     * @param x Positions of the atoms before computing the external forces,
     * with the following format:
     *
     * atom 1     atom 2     atom 3
     * x, y, z,   x, y, z,   x, y, z
     *
     * @param xprime Positions of the atoms after computing the external forces,
     * with the same format as x. When returning, it will contain the final
     * position of each atom.
     * @param vprime The atoms velocities. When returning, it will contain the
     * final velocity of each atom.
     * @param tol Maximum error allowed in each atom position.
     * @param maxit Maximum number of iterations the solver can execute.
     * @param deltat The time-step in picoseconds.
     * @param constraint_virial Update the virial?
     * @param virial sum r x m delta_r (virial)
     * @param constraint_velocities Update the velocities?
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param print_results Print accuracy results?
     * @return True if the constraints has been satisfied, false otherwise.
     */
    bool solve_newton(int time_step,
                      int &niters,
                      ArrayRef<const RVec> x,
                      ArrayRef<RVec> xprime,
                      ArrayRef<RVec> vprime,
                      real tol,
                      int maxit,
                      real deltat,
                      bool constraint_virial,
                      tensor virial,
                      bool constraint_velocities,
                      const t_pbc *pbc,
                      bool print_results);

    /**
     * Try to solve the bond constraints of MOL in at most MAXIT iterations of
     * the Simplified Newtons's method with a tolerance for each atom of TOL.
     *
     * @param time_step The time step.
     * @param niters When returning, NITERS will contain the total number of
     * iterations executed by the solver.
     *
     * @param x Positions of the atoms before computing the external forces,
     * with the following format:
     *
     * atom 1     atom 2     atom 3
     * x, y, z,   x, y, z,   x, y, z
     *
     * @param xprime Positions of the atoms after computing the external forces,
     * with the same format as x. When returning, it will contain the final
     * position of each atom.
     * @param vprime The atoms velocities. When returning, it will contain the
     * final velocity of each atom.
     * @param tol Maximum error allowed in each atom position.
     * @param maxit Maximum number of iterations the solver can execute.
     * @param deltat The time-step in picoseconds.
     * @param constraint_virial Update the virial?
     * @param virial sum r x m delta_r (virial)
     * @param constraint_velocities Update the velocities?
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param print_results Print accuracy results?
     * @return True if the constraints has been satisfied, false otherwise.
     */
    bool solve_simpl(int time_step,
                     int &niters,
                     ArrayRef<const RVec> x,
                     ArrayRef<RVec> xprime,
                     ArrayRef<RVec> vprime,
                     real tol,
                     int maxit,
                     real deltat,
                     bool constraint_virial,
                     tensor virial,
                     bool constraint_velocities,
                     const t_pbc *pbc,
                     bool print_results);

    /**
     * Try to solve the bond constraints of MOL in at most MAXIT iterations of
     * Forsgren's method with a tolerance for each atom of TOL.
     *
     * @param time_step The time step.
     * @param niters When returning, NITERS will contain the total number of
     * iterations executed by the solver.
     *
     * @param x Positions of the atoms before computing the external forces,
     * with the following format:
     *
     * atom 1     atom 2     atom 3
     * x, y, z,   x, y, z,   x, y, z
     *
     * @param xprime Positions of the atoms after computing the external forces,
     * with the same format as x. When returning, it will contain the final
     * position of each atom.
     * @param vprime The atoms velocities. When returning, it will contain the
     * final velocity of each atom.
     * @param tol Maximum error allowed in each atom position.
     * @param maxit Maximum number of iterations the solver can execute.
     * @param deltat The time-step in picoseconds.
     * @param constraint_virial Update the virial?
     * @param virial sum r x m delta_r (virial)
     * @param constraint_velocities Update the velocities?
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param print_results Print accuracy results?
     * @return True if the constraints has been satisfied, false otherwise.
     */
    bool solve_forsgren(int time_step,
                        int &niters,
                        ArrayRef<const RVec> x,
                        ArrayRef<RVec> xprime,
                        ArrayRef<RVec> vprime,
                        real tol,
                        int maxit,
                        real deltat,
                        bool constraint_virial,
                        tensor virial,
                        bool constraint_velocities,
                        const t_pbc *pbc,
                        bool print_results);


    /**
     * Try to solve the bond constraints of MOL in at most MAXIT iterations of
     * Quasi-Newton's methods with a tolerance for each atom of TOL.
     * This function is used to generate the data described in section 6.2.2.
     *
     * @param time_step The time step.
     * @param niters When returning, NITERS will contain the total number of
     * iterations executed by the solver.
     *
     * @param x Positions of the atoms before computing the external forces,
     * with the following format:
     *
     * atom 1     atom 2     atom 3
     * x, y, z,   x, y, z,   x, y, z
     *
     * @param xprime Positions of the atoms after computing the external forces,
     * with the same format as x. When returning, it will contain the final
     * position of each atom.
     * @param vprime The atoms velocities. When returning, it will contain the
     * final velocity of each atom.
     * @param tol Maximum error allowed in each atom position.
     * @param maxit Maximum number of iterations the solver can execute.
     * @param deltat The time-step in picoseconds.
     * @param constraint_virial Update the virial?
     * @param virial sum r x m delta_r (virial)
     * @param constraint_velocities Update the velocities?
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param print_results Print accuracy results?
     * @return True if the constraints has been satisfied, false otherwise.
     */
    bool solve_6_2_2(ArrayRef<const RVec> x,
                     ArrayRef<const RVec> xprime,
                     ArrayRef<const RVec> vprime,
                     real tol,
                     int maxit,
                     const tensor virial,
                     const t_pbc *pbc);

    /**
     * Try to solve the bond constraints of MOL in at most MAXIT iterations of
     * Newton's, Simplified-Newton's and Forsgren's (all-three) methods with a
     * tolerance for each atom of TOL. This function is used to generate the data
     * described in section 6.2.3.
     *
     * @param time_step The time step.
     * @param niters When returning, NITERS will contain the total number of
     * iterations executed by the solver.
     *
     * @param x Positions of the atoms before computing the external forces,
     * with the following format:
     *
     * atom 1     atom 2     atom 3
     * x, y, z,   x, y, z,   x, y, z
     *
     * @param xprime Positions of the atoms after computing the external forces,
     * with the same format as x. When returning, it will contain the final
     * position of each atom.
     * @param vprime The atoms velocities. When returning, it will contain the
     * final velocity of each atom.
     * @param tol Maximum error allowed in each atom position.
     * @param maxit Maximum number of iterations the solver can execute.
     * @param deltat The time-step in picoseconds.
     * @param constraint_virial Update the virial?
     * @param virial sum r x m delta_r (virial)
     * @param constraint_velocities Update the velocities?
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param print_results Print accuracy results?
     * @return True if the constraints has been satisfied, false otherwise.
     */
    bool solve_6_2_3(int time_step,
                     int &niters,
                     ArrayRef<const RVec> x,
                     ArrayRef<RVec> xprime,
                     ArrayRef<RVec> vprime,
                     real tol,
                     int maxit,
                     real deltat,
                     bool constraint_virial,
                     tensor virial,
                     bool constraint_velocities,
                     const t_pbc *pbc,
                     bool print_results);

private:
    molecule_t *mol;   // The molecule structure.

    // CRS sparse matrix.
    struct CRS {
        size_t n;   // rows and cols of the matrix (square matrix).

        // Indexes[i] is the index of the first entry in data, cols and rows of
        // the ith row.
        // Indexes[i + 1] is the last index + 1 of the first entry in data, cols
        // and rows of the ith row. Indexes contains n + 1 entries.
        std::vector<size_t> indexes;

        std::vector<size_t> diagonals;   // Indexes of the diagonal entries of
                                         // the matrix.

        std::vector<size_t> rows;   // The row of each index.
        std::vector<size_t> cols;   // The col of each index.

        std::vector<real> data;   // The nonzero entries of the matrix (also
                                  // fillins).

        std::vector<real> temp_work;   // A vector of n elements, used to store
                                       // temporal information.

        std::vector<real> weights;   // The weight of each nonzero entries of
                                     // the matrix (also fillins).
    };

    // START Experiments 6.2.2
    CRS A_newton;
    CRS A_quasi;

    std::ofstream reshis_file;
    std::ofstream normEk_file;
    std::ofstream rk_file;
    // END

    // START Experiments 6.2.3
    // Output files.
    std::ofstream newton_file;
    std::ofstream simpl_file;
    std::ofstream forsgren_file;

    // Times.
    int micros_newton;
    int micros_simpl;
    int micros_forsgren;
    // Number of iteration.
    int iters_newton;
    int iters_simpl;
    int iters_forsgren;
    std::ofstream times_file;

    // Coefficient matrices.
    CRS A_simpl;
    CRS A_forsgren_init;
    CRS A_forsgren;

    // Temporal arrays for Forsgren's method.
    std::vector<real> dk;
    std::vector<real> fk;
    std::vector<real> Fpzkdk;
    // END

    // Independent term vector. This vector will hold the lagrange multipliers
    // at the end of an ILVES iteration.
    std::vector<real> g;

    // The current approximation of the Lagrange multipliers.
    std::vector<real> current_lagr;
    // The final approximation of the Lagrange multipliers for the newton
    // method.
    std::vector<real> lagr_newton;

    // x_ab[i][XX/YY/ZZ] contains the vector from atom a to b (using x),
    // the atoms that are part of the ith bond.
    std::vector<RVec> x_ab;
    // x_ab[i][XX/YY/ZZ] contains the vector from atom a to b (using xprime),
    // the atoms that are part of the ith bond.
    std::vector<RVec> xprime_ab;

    /**
     * Apply a reordering to the bonds of MOL to reduce the number of fillins
     * generated when applying an LU/Cholesky factorization.
     *
     */
    void minimal_degree_graph_reordering();

    /**
     * Copy xprime, vprime and virial to xprime_back, vprime_back and
     * virial_back.
     *
     * @param xprime Array of atomic positions.
     * @param vprime Array of atomic velocities.
     * @param virial Virial tensor.
     * @param xprime_back Copy of xprime.
     * @param vprime_back Copy of vprime.
     * @param virial_back Copy of virial.
     */
    void backup(ArrayRef<const RVec> xprime,
                ArrayRef<const RVec> vprime,
                const tensor virial,
                ArrayRef<RVec> xprime_back,
                ArrayRef<RVec> vprime_back,
                tensor virial_back);

    /**
     * Print accuray results to file.
     *
     * @param file The output file.
     * @param time_step The time-step.
     * @param iter The iteration number.
     * @param tau The maximum relative bond length violation.
     */
    void print_to_file(std::ofstream &file, int time_step, int iter, real tau);

    /**
     * Subtraction of two rvec taking the PBC into account if pbc is not null.
     *
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param vin1 First source rvec.
     * @param vin2 Second source rvec.
     * @param vout Destination rvec.
     */
    void pbc_rvec_sub(const t_pbc *pbc, const rvec vin1, const rvec vin2, rvec vout);

    /**
     * Constructs g (the right-hand side of the linear system).
     * Also updates xprime_ab and x_ab (if compute_x_ab).
     * Returns the largest relative (square) bond length violation.
     *
     * @param pbc The PBC (container) information. Null if there is no PBC
     * information.
     * @param x Positions before computing the external forces.
     * @param xprime Positions of the atoms after computing the external forces
     * (current positions of the atoms).
     * @param compute_x_ab If true, x_ab will be updated.
     * @return real The largest relative (square) bond length violation.
     */
    real make_g(const t_pbc *pbc,
                ArrayRef<const RVec> x,
                ArrayRef<const RVec> xprime,
                bool compute_x_ab);

    /**
     * Construct a the left-hand-size of the linear system given a sparse
     * matrix.
     *
     * @param A The sparse matrix.
     */
    void make_A_asym(CRS &A);

    /**
     * Construct a the left-hand-size of the linear system given an
     * upper-triangular sparse matrix.
     *
     * @param A The upper-triangular sparse matrix.
     */
    void make_A_sym(CRS &A);

    /**
     * Performs the sparse LU decomposition of the sparse matrix A. L and U will
     * overwrite A.
     *
     * @param A The sparse matrix.
     *
     */
    void LU(CRS &A);

    /**
     * Performs the sparse Cholesky decomposition of the upper-triangular sparse
     * matrix A. L and LT will overwrite A.
     *
     * @param A The upper-triangular sparse matrix.
     *
     */
    void cholesky(CRS &A);

    /**
     * Solve Lx = g using forward substitution, overwriting g with the solution
     * x.
     *
     * @param LU The sparse matrix.
     *
     */
    void forward_LU(const CRS &LU);

    /**
     * Solve Uy = g using backward substitution, overwriting g with the
     * solution y.
     *
     * @param LU The sparse matrix.
     *
     */
    void backward_LU(const CRS &LU);

    /**
     * Solve Lx = g using forward substitution, overwriting g with the solution
     * x.
     *
     * @param LLT The upper-triangular sparse matrix.
     *
     */
    void forward_cholesky(const CRS &LLT);

    /**
     * Solve LTy = g using backward substitution, overwriting g with the
     * solution y.
     *
     * @param LLT The upper-triangular sparse matrix.
     *
     */
    void backward_cholesky(const CRS &LLT);

    /**
     * For each bond b, that joins atoms a and b:
     *     XPRIME[a] += lagr_correction[b] * (X>[b] - X>[a]) * (1 / mass[a]);
     *     XPRIME[b] -= lagr_correction[b] * (X>[b] - X>[a]) * (1 / mass[a]);
     *
     * @param xprime Positions of the atoms after computing the external forces
     * (current positions of the atoms).
     */
    void update_positions(ArrayRef<RVec> xprime);

    /**
     * For each bond b that joins atoms a and b:
     *
     * for d1 in DIM
     *     for d2 in DIM
     *         VIRIAL[d1][d2] -= -lagr[b] * (X>[b] - X>[a])[d1] * (X>[b]
     * - X>[a])[d2]
     *
     * @param virial The virial (a DIMxDIM matrix).
     */
    void update_virial(tensor virial);

    /**
     * For each bond b that joins atoms a and b:
     *     V[a] += current_lagr[b] * (X>[b] - X>[a]) * (1 / mass[a]) * INVDT;
     *     V[b] -= current_lagr[b] * (X>[b] - X>[a]) * (1 / mass[a]) * INVDT;
     *
     * @param vprime The atoms velocities.
     * @param invdt Inverse of the time-step in picoseconds.
     */
    void update_velocities(ArrayRef<RVec> vprime, real invdt);
};

}   // namespace gmx
