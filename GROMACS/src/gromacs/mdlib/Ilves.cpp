#include "Ilves.h"

#include <stdio.h>

#include <algorithm>
#include <chrono>
#include <iomanip>

#include "gromacs/math/vec.h"

static real norm2(const gmx::ArrayRef<const gmx::RVec> &v1,
                  const gmx::ArrayRef<const gmx::RVec> &v2) {
    real norm = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        for (size_t d = 0; d < DIM; ++d) {
            norm += std::pow(v1[i][d] - v2[i][d], 2);
        }
    }
    norm = std::sqrt(norm);

    return norm;
}

static real norm2(const gmx::ArrayRef<const gmx::RVec> &v) {
    real norm = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t d = 0; d < DIM; ++d) {
            norm += std::pow(v[i][d], 2);
        }
    }
    norm = std::sqrt(norm);

    return norm;
}

static real norm2(const std::vector<real> &v1, const std::vector<real> &v2) {
    real norm = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        norm += std::pow(v1[i] - v2[i], 2);
    }
    norm = std::sqrt(norm);

    return norm;
}

static real norm2(const std::vector<real> &v) {
    real norm = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        norm += std::pow(v[i], 2);
    }
    norm = std::sqrt(norm);

    return norm;
}

gmx::Ilves::Ilves(molecule_t *const mol) : mol(mol) {
    minimal_degree_graph_reordering();

    x_ab.resize(mol->n);
    xprime_ab.resize(mol->n);

    current_lagr.resize(mol->n);
    lagr_newton.resize(mol->n);

    g.resize(mol->n);

    A_newton.n = mol->n;

    A_newton.temp_work.resize(A_newton.n, 0);

    // Construct the CRS representation of A.
    // Asymmetric matrix for Newton's and Forsgren's method.

    size_t idx = 0;

    A_newton.indexes.push_back(idx);

    for (size_t row = 0; row < A_newton.n; ++row) {
        // Loop over the nonzeros in the rowth row
        for (int l = mol->bond_graph->xadj[row];
             l < mol->bond_graph->xadj[row + 1];
             ++l, ++idx) {
            // Determine bo of the idxth entry.
            const auto col = mol->bond_graph->adj[l];

            A_newton.cols.push_back(col);
            A_newton.rows.push_back(row);

            A_newton.weights.push_back(mol->weights[idx]);

            if (col == row) {
                A_newton.diagonals.push_back(idx);
            }
        }
        A_newton.indexes.push_back(idx);
    }

    A_newton.data.resize(idx);

    A_forsgren_init = A_newton;
    A_forsgren = A_newton;

    // Symmetric matrix for the simpl-Newton method.

    A_simpl.n = mol->n;

    A_simpl.temp_work.resize(A_simpl.n, 0);

    // Construct the Upper triangular CRS representation of AS.

    idx = 0;
    size_t weights_idx = 0;

    A_simpl.indexes.push_back(idx);

    for (size_t row = 0; row < A_simpl.n; ++row) {
        // Loop over the nonzeros in the rowth row
        for (int l = mol->bond_graph->xadj[row];
             l < mol->bond_graph->xadj[row + 1];
             ++l, ++weights_idx) {
            // Determine the col of the lth entry.
            const auto col = mol->bond_graph->adj[l];

            // We only need the sparsity pattern of the upper triangular part
            // of the adjacency matrix.
            if (col >= row) {
                A_simpl.rows.push_back(row);
                A_simpl.cols.push_back(col);

                A_simpl.weights.push_back(mol->weights[weights_idx]);

                ++idx;
            }
        }
        A_simpl.indexes.push_back(idx);
    }

    A_simpl.data.resize(idx);

    A_quasi = A_simpl;

    // Temporal arrays for Forsgren's method.
    dk.resize(g.size());
    fk.resize(g.size());
    Fpzkdk.resize(g.size());

    reshis_file.open("rehis.txt");
    normEk_file.open("normEk.txt");
    rk_file.open("rk.txt");

    newton_file.open("newton.txt");
    simpl_file.open("simpl_newton.txt");
    forsgren_file.open("forsgren.txt");

    newton_file << "time-step\titeration\t"
                   "||lagr-newton||\ttau\t"
                   "||lagr-correction-iter||\t||current-lagr||\t"
                   "||lagr-newton - current_lagr||\t"
                   "||lagr-newton - current-lagr||/||lagr-newton||";

    simpl_file << "time-step\titeration\t"
                  "||lagr-newton||\ttau\t"
                  "||lagr-correction-iter||\t||current-lagr||\t"
                  "||lagr-newton - current_lagr||\t"
                  "||lagr-newton - current-lagr||/||lagr-newton||";

    forsgren_file << "time-step\titeration\t"
                     "||lagr-newton||\ttau\t"
                     "||lagr-correction-iter||\t||current-lagr||\t"
                     "||lagr-newton - current_lagr||\t"
                     "||lagr-newton - current-lagr||/||lagr-newton||";

    times_file.open("times.txt");

    micros_newton = 0;
    micros_simpl = 0;
    micros_forsgren = 0;

    iters_newton = 0;
    iters_simpl = 0;
    iters_forsgren = 0;
}

gmx::Ilves::~Ilves() {
    times_file << "newton\tsimpl\tforsgren\n";
    times_file << micros_newton << "\t" << micros_simpl << "\t"
               << micros_forsgren << "\n";
    times_file << iters_newton << "\t" << iters_simpl << "\t" << iters_forsgren
               << "\n";
}

void gmx::Ilves::print_to_file(std::ofstream &file,
                               const int time_step,
                               const int iter,
                               const real tau) {
    if (iter < 0) {
        file << time_step << "\t"
             << "initial"
             << "\t"
             << "-"
             << "\t" << tau << "\t"
             << "-"
             << "\t"
             << "-"
             << "\t"
             << "-"
             << "\t"
             << "-"
             << "\n";
    }
    else {
        real norm_lagr_newton = ::norm2(lagr_newton);
        real norm_lagr_correction = ::norm2(g);
        real norm_current_lagr = ::norm2(current_lagr);
        real norm_lagr_newton_current_lagr = ::norm2(lagr_newton, current_lagr);

        file << time_step << "\t" << iter << "\t" << norm_lagr_newton << "\t"
             << tau << "\t" << norm_lagr_correction << "\t" << norm_current_lagr
             << "\t" << norm_lagr_newton_current_lagr << "\t"
             << norm_lagr_newton_current_lagr / norm_lagr_newton << "\n";
    }
}

void gmx::Ilves::minimal_degree_graph_reordering() {
    /*
     * Reorder the nodes of the graph to reduce the fillins.
     */
    int *permutation = (int *)malloc(mol->n * sizeof(permutation[0]));

    graph_minimal_degree(mol->bond_graph, permutation);
    renumber_bonds(mol, permutation);

    free(permutation);

    /*
     * Compute the fillin graph.
     */

    // Allocate space for the new bond graphs
    graph_t *full_fill_graph = (graph_t *)malloc(sizeof(*full_fill_graph));

    // Compute the fill.
    graph_compute_fill(mol->bond_graph, full_fill_graph);

    graph_free(mol->bond_graph);
    free(mol->bond_graph);

    mol->bond_graph = full_fill_graph;

    make_weights(mol, mol->bond_graph);
}

void gmx::Ilves::pbc_rvec_sub(const t_pbc *const pbc,
                              const rvec vin1,
                              const rvec vin2,
                              rvec vout) {
    if (pbc) {
        pbc_dx_aiuc(pbc, vin1, vin2, vout);
    }
    else {
        rvec_sub(vin1, vin2, vout);
    }
}

bool gmx::Ilves::solve_6_2_2(const ArrayRef<const RVec> x,
                             const ArrayRef<const RVec> xprime,
                             const ArrayRef<const RVec> vprime,
                             const real tol,
                             const int maxit,
                             const tensor virial,
                             const t_pbc *const pbc) {
    std::vector<real> z = current_lagr;

    std::vector<RVec> xprime_initial(mol->m);
    std::vector<RVec> vprime_initial(mol->m);
    tensor virial_initial;

    backup(xprime, vprime, virial, xprime_initial, vprime_initial, virial_initial);

    std::vector<RVec> xprime_back(mol->m);
    std::vector<RVec> vprime_back(mol->m);
    tensor virial_back;

    // Get the reference solution.
    bool success;
    {
        int niters;
        backup(xprime, vprime, virial, xprime_back, vprime_back, virial_back);
        success = solve_newton(0,
                               niters,
                               x,
                               xprime_back,
                               vprime_back,
                               tol,
                               maxit,
                               1,
                               false,
                               virial_back,
                               false,
                               pbc,
                               false);
        std::swap(z, current_lagr);
    }

    // Compute g(x). Update x_ab and xprime_ab.
    real tau = make_g(pbc, x, xprime_initial, true);

    real max_g = 0;
    for (size_t row = 0; row < g.size(); ++row) {
        max_g = std::max(max_g, std::abs(g[row]));
    }

    // Corrections.
    std::vector<real> sk(mol->n);   // Newton method
    std::vector<real> tk(mol->n);   // Quasi-newton method

    // Lagrange multipliers.
    std::vector<real> xk(mol->n, 0);   // Newton method
    std::vector<real> yk(mol->n, 0);   // Quasi-newton method

    real norm_z = ::norm2(z);

    fprintf(stderr, "max-rel-error-squared-quasi-init\n");
    fprintf(stderr, "%E\n", tau);
    fprintf(stderr, "||z||\n");
    fprintf(stderr, "%E\n", norm_z);
    fprintf(stderr,
            "k\t||xk||\t||sk||\t||yk||\t||tk||\t||Ek||\t||z - yk||\t||z - "
            "yk||/||z||\tmax-rel-error-squared-quasi-k\n");

    // Do at most MAXIT Newton steps.
    for (int i = 0; i < maxit && tol < tau; ++i) {
        /*
         * Newton method.
         */

        // Construct A.
        make_A_asym(A_newton);

        // Solve the linear system.
        LU(A_newton);
        forward_LU(A_newton);
        backward_LU(A_newton);

        // Store the correction of the lagrange multiplier
        sk = g;

        // Restore gx.
        make_g(pbc, x, xprime_initial, false);

        /*
         * Quasi-newton method.
         */

        // Construct A.
        make_A_sym(A_quasi);

        // Solve the linear system.
        cholesky(A_quasi);
        forward_cholesky(A_quasi);
        backward_cholesky(A_quasi);

        // Store the correction of the lagrange multiplier
        tk = g;

        for (size_t row = 0; row < g.size(); ++row) {
            current_lagr[row] = (i == 0) ? g[row] : current_lagr[row] + g[row];
        }

        // Print results
        real norm_xk = ::norm2(xk);
        real norm_sk = ::norm2(sk);

        real norm_yk = ::norm2(yk);
        real norm_tk = ::norm2(tk);

        real norm_sk_tk = ::norm2(sk, tk);

        real norm_z_yk = ::norm2(z, yk);

        real Ek = norm_sk_tk / norm_sk;

        reshis_file << std::setprecision(17) << max_g << "\t";
        normEk_file << std::setprecision(17) << Ek << "\t";
        rk_file << std::setprecision(17) << norm_z_yk / norm_z << "\t";

        fprintf(stderr,
                "%d\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n",
                i,
                norm_xk,
                norm_sk,
                norm_yk,
                norm_tk,
                Ek,
                norm_z_yk,
                norm_z_yk / norm_z,
                tau);

        // Update the lagrange multipliers.
        for (size_t bond = 0; bond < g.size(); ++bond) {
            yk[bond] += tk[bond];
            xk[bond] += sk[bond];
        }

        // Update the atoms position.
        update_positions(xprime_initial);

        // Compute g(x). Update xprime_ab.
        tau = make_g(pbc, x, xprime_initial, false);

        max_g = 0;
        for (size_t row = 0; row < g.size(); ++row) {
            max_g = std::max(max_g, std::abs(g[row]));
        }
    }
    // error = tau;

    real norm_z_yk = ::norm2(z, yk);

    reshis_file << std::setprecision(17) << max_g << "\t";
    rk_file << std::setprecision(17) << norm_z_yk / norm_z << "\t";

    reshis_file << "\n";
    normEk_file << "\n";
    rk_file << "\n";

    reshis_file.flush();
    normEk_file.flush();
    rk_file.flush();

    return tau < tol;
}

bool gmx::Ilves::solve_newton(const int time_step,
                              int &niters,
                              const ArrayRef<const RVec> x,
                              const ArrayRef<RVec> xprime,
                              const ArrayRef<RVec> vprime,
                              const real tol,
                              const int maxit,
                              const real deltat,
                              const bool constraint_virial,
                              tensor virial,
                              const bool constraint_velocities,
                              const t_pbc *const pbc,
                              const bool print_results) {
    niters = 0;

    // Compute g(x). Update x_ab and xprime_ab.
    real tau = make_g(pbc, x, xprime, true);
    if (print_results) {
        print_to_file(newton_file, time_step, -1, tau);
    }

    // Do at most MAXIT Newton steps.
    for (int i = 0; i < maxit && tol < tau; ++i) {
        ++niters;

        // Construct A.
        make_A_asym(A_newton);

        // Solve the linear system.
        LU(A_newton);
        forward_LU(A_newton);
        backward_LU(A_newton);

        // Update the lagrange multipliers.
        for (size_t bond = 0; bond < g.size(); ++bond) {
            current_lagr[bond] = (i == 0) ? g[bond] : current_lagr[bond] + g[bond];
        }

        update_positions(xprime);

        // Compute g(x). Update xprime_ab.
        tau = make_g(pbc, x, xprime, false);

        if (print_results) {
            print_to_file(newton_file, time_step, i, tau);
        }
    }

    // Update the virial.
    if (constraint_virial) {
        update_virial(virial);
    }

    // Update the velocities.
    if (constraint_velocities) {
        update_velocities(vprime, 1 / deltat);
    }

    return tau < tol;
}

bool gmx::Ilves::solve_simpl(const int time_step,
                             int &niters,
                             const ArrayRef<const RVec> x,
                             const ArrayRef<RVec> xprime,
                             const ArrayRef<RVec> vprime,
                             const real tol,
                             const int maxit,
                             const real deltat,
                             const bool constraint_virial,
                             tensor virial,
                             const bool constraint_velocities,
                             const t_pbc *const pbc,
                             const bool print_results) {
    niters = 0;

    // Compute g(x). Update x_ab and xprime_ab.
    real tau = make_g(pbc, x, xprime, true);
    if (print_results) {
        print_to_file(simpl_file, time_step, -1, tau);
    }

    // Construct A.
    make_A_sym(A_simpl);
    // Solve the linear system.
    cholesky(A_simpl);

    // Do at most MAXIT Newton steps.
    for (int i = 0; i < maxit && tol < tau; ++i) {
        ++niters;

        forward_cholesky(A_simpl);
        backward_cholesky(A_simpl);

        // Update the lagrange multipliers.
        for (size_t bond = 0; bond < g.size(); ++bond) {
            current_lagr[bond] = (i == 0) ? g[bond] : current_lagr[bond] + g[bond];
        }

        update_positions(xprime);

        // Compute g(x). Update xprime_ab.
        tau = make_g(pbc, x, xprime, false);

        if (print_results) {
            print_to_file(simpl_file, time_step, i, tau);
        }
    }

    // Update the virial.
    if (constraint_virial) {
        update_virial(virial);
    }

    // Update the velocities.
    if (constraint_velocities) {
        update_velocities(vprime, 1 / deltat);
    }

    return tau < tol;
}

bool gmx::Ilves::solve_forsgren(const int time_step,
                                int &niters,
                                const ArrayRef<const RVec> x,
                                const ArrayRef<RVec> xprime,
                                const ArrayRef<RVec> vprime,
                                const real tol,
                                const int maxit,
                                const real deltat,
                                const bool constraint_virial,
                                tensor virial,
                                const bool constraint_velocities,
                                const t_pbc *const pbc,
                                const bool print_results) {
    niters = 0;
    // Compute g(x). Update x_ab and xprime_ab.
    real tau = make_g(pbc, x, xprime, true);
    if (print_results) {
        print_to_file(forsgren_file, time_step, -1, tau);
    }

    make_A_asym(A_forsgren_init);   // F'(z0) for Forsgren method.
    LU(A_forsgren_init);

    // Do at most MAXIT Newton steps.
    for (int i = 0; i < maxit && tol < tau; ++i) {
        ++niters;

        int mk = std::pow(2, i);

        make_A_asym(A_forsgren);   // A = F'(zk)

        std::swap(fk, g);   // F(zk)

        for (int j = 0; j < mk; ++j) {
            if (j > 0) {
                std::fill(Fpzkdk.begin(), Fpzkdk.end(), 0);

                for (size_t idx = 0; idx < A_forsgren.indexes.back(); ++idx) {
                    const auto row = A_forsgren.rows[idx];
                    const auto col = A_forsgren.cols[idx];

                    Fpzkdk[row] += A_forsgren.data[idx] * dk[col];
                }
            }

            for (size_t bond = 0; bond < g.size(); ++bond) {
                g[bond] = fk[bond] - ((j == 0) ? 0 : Fpzkdk[bond]);
            }

            forward_LU(A_forsgren_init);
            backward_LU(A_forsgren_init);

            for (size_t bond = 0; bond < g.size(); ++bond) {
                dk[bond] = (j == 0) ? g[bond] : dk[bond] + g[bond];
            }
        }

        std::swap(g, dk);

        // Update the lagrange multipliers.
        for (size_t bond = 0; bond < g.size(); ++bond) {
            current_lagr[bond] = (i == 0) ? g[bond] : current_lagr[bond] + g[bond];
        }

        update_positions(xprime);

        // Compute g(x). Update xprime_ab.
        tau = make_g(pbc, x, xprime, false);

        if (print_results) {
            print_to_file(forsgren_file, time_step, i, tau);
        }
    }

    // Update the virial.
    if (constraint_virial) {
        update_virial(virial);
    }

    // Update the velocities.
    if (constraint_velocities) {
        update_velocities(vprime, 1 / deltat);
    }

    return tau < tol;
}

bool gmx::Ilves::solve_6_2_3(int time_step,
                             int &,   // niters, unused.
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
                             bool print_results) {
    if (print_results) {
        newton_file << "\n";
        simpl_file << "\n";
        forsgren_file << "\n";
    }

    std::vector<RVec> xprime_initial(mol->m);
    std::vector<RVec> vprime_initial(mol->m);
    tensor virial_initial;

    backup(xprime, vprime, virial, xprime_initial, vprime_initial, virial_initial);

    std::vector<RVec> xprime_back(mol->m);
    std::vector<RVec> vprime_back(mol->m);
    tensor virial_back;

    //--------------------------------------------------------------------------
    // Get the reference solution.
    bool success;
    {
        int niters;
        backup(xprime, vprime, virial, xprime_back, vprime_back, virial_back);
        success = solve_newton(time_step,
                               niters,
                               x,
                               xprime,
                               vprime,
                               tol,
                               maxit,
                               deltat,
                               constraint_virial,
                               virial,
                               constraint_velocities,
                               pbc,
                               false);
        std::swap(lagr_newton, current_lagr);
    }

    // Solve with all methods, and compare to Newton's solution.
    //--------------------------------------------------------------------------
    {
        backup(xprime_initial,
               vprime_initial,
               virial_initial,
               xprime_back,
               vprime_back,
               virial_back);

        int niters;
        auto start = std::chrono::steady_clock::now();
        solve_newton(time_step,
                     niters,
                     x,
                     xprime_back,
                     vprime_back,
                     tol,
                     maxit,
                     deltat,
                     constraint_virial,
                     virial_back,
                     constraint_velocities,
                     pbc,
                     print_results);
        auto end = std::chrono::steady_clock::now();
        micros_newton += std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                             .count();
        iters_newton += niters;
    }

    //--------------------------------------------------------------------------
    {
        backup(xprime_initial,
               vprime_initial,
               virial_initial,
               xprime_back,
               vprime_back,
               virial_back);

        int niters;
        auto start = std::chrono::steady_clock::now();
        solve_simpl(time_step,
                    niters,
                    x,
                    xprime_back,
                    vprime_back,
                    tol,
                    maxit,
                    deltat,
                    constraint_virial,
                    virial_back,
                    constraint_velocities,
                    pbc,
                    print_results);
        auto end = std::chrono::steady_clock::now();
        micros_simpl += std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                            .count();
        iters_simpl += niters;
    }
    //--------------------------------------------------------------------------
    {
        backup(xprime_initial,
               vprime_initial,
               virial_initial,
               xprime_back,
               vprime_back,
               virial_back);

        int niters;
        auto start = std::chrono::steady_clock::now();
        solve_forsgren(time_step,
                       niters,
                       x,
                       xprime_back,
                       vprime_back,
                       tol,
                       maxit,
                       deltat,
                       constraint_virial,
                       virial_back,
                       constraint_velocities,
                       pbc,
                       print_results);
        auto end = std::chrono::steady_clock::now();
        micros_forsgren += std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                               .count();
        iters_forsgren += niters;
    }

    //--------------------------------------------------------------------------

    // if (print_results) {
    //     newton_file.flush();
    //     simpl_file.flush();
    //     forsgren_file.flush();
    // }

    return success;
}

real gmx::Ilves::make_g(const t_pbc *const pbc,
                        const ArrayRef<const RVec> x,
                        const ArrayRef<const RVec> xprime,
                        const bool compute_x_ab) {
    // Initialize the largest relative error
    real rel = 0;

    // Loop over the n constraints
    for (size_t row = 0; row < g.size(); ++row) {
        // Isolate the atoms which partake in bond row
        const auto a = mol->bonds[2 * row];
        const auto b = mol->bonds[2 * row + 1];

        // Compute the vectors from atom a to atom b.
        if (compute_x_ab) {
            pbc_rvec_sub(pbc, x[b], x[a], x_ab[row].as_vec());
        }
        pbc_rvec_sub(pbc, xprime[b], xprime[a], xprime_ab[row].as_vec());

        // Compute the square of the length of r
        const auto scalar = iprod(xprime_ab[row].as_vec(), xprime_ab[row].as_vec());

        // Compute the constraint violation
        g[row] = 0.5 * (scalar - mol->sigma2[row]);

        // Update the relative error
        rel = std::max(rel, std::abs(g[row]) / mol->sigma2[row]);
    }

    // Return the largest relative (square) bond length violation
    return rel;
}

void gmx::Ilves::make_A_asym(CRS &A) {
    for (size_t idx = 0; idx < A.indexes.back(); ++idx) {
        const auto row = A.rows[idx];
        const auto col = A.cols[idx];

        // Compute the inner product scalar.
        const auto scalar = iprod(xprime_ab[row].as_vec(), x_ab[col].as_vec());

        A.data[idx] = scalar * A.weights[idx];
    }
}

void gmx::Ilves::make_A_sym(CRS &A) {
    for (size_t idx = 0; idx < A.indexes.back(); ++idx) {
        const auto row = A.rows[idx];
        const auto col = A.cols[idx];

        // Compute the inner product scalar.
        const auto scalar = iprod(x_ab[row].as_vec(), x_ab[col].as_vec());

        A.data[idx] = scalar * A.weights[idx];
    }
}

void gmx::Ilves::LU(CRS &A) {
    // Loop over the first n-1 columns of the matrix.
    for (size_t col = 0; col < A.n - 1; ++col) {
        // Isolate the diagonal entry A[col,col].
        const auto pivot = A.data[A.diagonals[col]];

        // Process all *relevant* rows below row col
        for (size_t idx = A.diagonals[col] + 1; idx < A.indexes[col + 1]; ++idx) {
            // WARNING:
            // This is correct only when the matrix has symmetric structure.

            // Isolate the relevant row for readability.
            const auto row = A.cols[idx];

            // Expand A(row,:) into temp_work.
            for (size_t s = A.indexes[row]; s < A.indexes[row + 1]; ++s) {
                A.temp_work[A.cols[s]] = A.data[s];
            }

            /*
             * Perform the Gaussian elimination.
             */

            A.temp_work[col] /= pivot;

            const auto f = A.temp_work[col];

            // Update the row (A[row, t] -= f * A[col, t]).
            for (size_t t = A.diagonals[col] + 1; t < A.indexes[col + 1]; ++t) {
                A.temp_work[A.cols[t]] -= f * A.data[t];
            }

            // Compress temp_work into A[row,:]. Also reset temp_work to avoid
            // pollution in the following iterations.
            for (size_t s = A.indexes[row]; s < A.indexes[row + 1]; ++s) {
                A.data[s] = A.temp_work[A.cols[s]];

                A.temp_work[A.cols[s]] = 0;
            }
        }
    }
}

void gmx::Ilves::cholesky(CRS &A) {
    // Loop over the first n rows of the matrix.
    for (size_t row = 0; row < A.n; ++row) {
        // Take the square root of the diagonal entry,
        // A(row,row):=sqrt(A(row,row))
        const auto diag = std::sqrt(A.data[A.indexes[row]]);
        A.data[A.indexes[row]] = diag;

        /* In the next loop we scale the upperdiagonal entries of the rowth row
           with A(row,row)
              A(row,row+1:n) /= A(row,row).
           Moreover, we expand the compressed representation of the updated
           entries into the array WORK
        */
        for (size_t s = A.indexes[row] + 1; s < A.indexes[row + 1]; ++s) {
            // Scale the subdiagonal entry
            A.data[s] /= diag;
            // Copy it into the auxillary vector work
            A.temp_work[A.cols[s]] = A.data[s];
        }

        /* At this point we have to perform a rank 1 update of the submatrix in
           the current upper right corner, specifically

           A(row+1:n,row+1:n) -= A(row+1:n,row) * A(row+1:n,j)'

           Only the lower triangular part matters. In MATLAB we would write

             for s=row+1:n
                alpha=A(row,s);
                A(s,s:n) -= alpha*A(row,s:n);
             end

          which stresses the fact that only nonzero entries alpha = A(row,s)
          matter.

        */

        // Loop over the nonzero upperdiagonal elements of the colth row.
        for (size_t s = A.indexes[row] + 1; s < A.indexes[row + 1]; ++s) {
            // Isolate the nonzero entry
            const auto alpha = A.data[s];
            // Isolate the index of the row which we have to update
            const auto col = A.cols[s];
            // Update the upperdiagonal part of this row
            for (size_t t = A.indexes[col]; t < A.indexes[col + 1]; ++t) {
                A.data[t] -= alpha * A.temp_work[A.cols[t]];
            }
        }

        /* At this point the jth rank 1 update is complete. We must clear the
           nonzero entries of WORK, so that it does not contain garbage during
           the next iteration.

           Notice that we do not bother to fill the entire array with zeros.
           We only kill those which could be nonzero.
        */
        for (size_t s = A.indexes[row] + 1; s < A.indexes[row + 1]; ++s) {
            A.temp_work[A.cols[s]] = 0;
        }
    }
}

void gmx::Ilves::forward_LU(const CRS &LU) {
    for (size_t row = 0; row < LU.n; ++row) {
        for (size_t idx = LU.indexes[row]; idx < LU.diagonals[row]; ++idx) {
            const auto f = LU.data[idx];

            g[row] -= f * g[LU.cols[idx]];
        }
    }
}

void gmx::Ilves::backward_LU(const CRS &LU) {
    for (size_t row = LU.n - 1; row != -1; --row) {
        for (size_t idx = LU.indexes[row + 1] - 1; idx > LU.diagonals[row]; --idx) {
            const auto f = LU.data[idx];

            g[row] -= f * g[LU.cols[idx]];
        }

        const auto pivot = LU.data[LU.diagonals[row]];

        g[row] /= pivot;
    }
}

void gmx::Ilves::forward_cholesky(const CRS &LLT) {
    // Loop over the first n rows of the matrix.
    for (size_t row = 0; row < LLT.n; ++row) {
        // Isolate the rowth diagonal entry of L
        const auto diag = LLT.data[LLT.indexes[row]];
        // Divide with diagonal entry
        g[row] /= diag;
        /* At this point x[row] has been computed and we must eliminate it from
           equations row+1, row+2, ..., n-1.
        */
        // Loop over the strictly upperdiagonal, nonzero entries of the rowth
        // row
        for (size_t s = LLT.indexes[row] + 1; s < LLT.indexes[row + 1]; ++s) {
            // Removing the influence of x[row] from row rows[diag] of the RHS
            g[LLT.cols[s]] -= LLT.data[s] * g[row];
        }
    }
}

void gmx::Ilves::backward_cholesky(const CRS &LLT) {
    // Loop over the first n rows of the matrix backwards.
    for (size_t row = LLT.n - 1; row != -1; --row) {
        // Isolate the rowth diagonal element
        const auto diag = LLT.data[LLT.indexes[row]];
        // Remove the contribution of the variables x[row+1], ..., x[n-1]
        for (size_t s = LLT.indexes[row] + 1; s < LLT.indexes[row + 1]; ++s) {
            g[row] -= LLT.data[s] * g[LLT.cols[s]];
        }
        // Divide with the diagonal element
        g[row] /= diag;
    }
}

void gmx::Ilves::update_positions(const ArrayRef<RVec> xprime) {
    for (int row = 0; row < mol->n; ++row) {
        const auto a = mol->bonds[2 * row];
        const auto b = mol->bonds[2 * row + 1];

        // Update the components of y which correspond to atom a
        for (size_t d = 0; d < DIM; ++d) {
            xprime[a][d] += g[row] * x_ab[row][d] * mol->invmass[a];
        }

        // Update the components of y which correspond to atom b
        for (size_t d = 0; d < DIM; ++d) {
            xprime[b][d] -= g[row] * x_ab[row][d] * mol->invmass[b];
        }
    }
}

void gmx::Ilves::update_virial(tensor virial) {
    for (int row = 0; row < mol->n; ++row) {
        for (size_t d1 = 0; d1 < DIM; ++d1) {
            const auto tmp = -current_lagr[row] * x_ab[row][d1];

            for (size_t d2 = 0; d2 < DIM; d2++) {
                virial[d1][d2] -= tmp * x_ab[row][d2];
            }
        }
    }
}

void gmx::Ilves::update_velocities(const ArrayRef<RVec> vprime, const real invdt) {
    for (int row = 0; row < mol->n; ++row) {
        auto a = mol->bonds[2 * row];
        auto b = mol->bonds[2 * row + 1];   // Atoms involved in the bond row

        for (size_t d = 0; d < DIM; ++d) {
            vprime[a][d] += current_lagr[row] * x_ab[row][d] * mol->invmass[a] *
                            invdt;
        }

        for (size_t d = 0; d < DIM; ++d) {
            vprime[b][d] -= current_lagr[row] * x_ab[row][d] * mol->invmass[b] *
                            invdt;
        }
    }
}

void gmx::Ilves::backup(const ArrayRef<const RVec> xprime,
                        const ArrayRef<const RVec> vprime,
                        const tensor virial,
                        ArrayRef<RVec> xprime_back,
                        ArrayRef<RVec> vprime_back,
                        tensor virial_back) {
    for (size_t atom = 0; atom < mol->m; ++atom) {
        for (int d = 0; d < DIM; ++d) {
            xprime_back[atom][d] = xprime[atom][d];
            vprime_back[atom][d] = vprime[atom][d];
        }
    }

    for (int d1 = 0; d1 < DIM; ++d1) {
        for (int d2 = 0; d2 < DIM; ++d2) {
            virial_back[d1][d2] = virial[d1][d2];
        }
    }
}
