# [PAPER] Newton's method revisited: How accurate do we have to be?

This repository contains the code used to generate the data for Figures 1, 2, 3, 4, 5, and 6 of the paper titled [Newton's method revisited: How accurate do we have to be?](TODO). What can you find in this repository?

1. The MATLAB codes necessary to compute square roots and generate Figure 1. MATLAB is a commercial product and the property of MathWorks, hence it cannot be included in this repository. MATLAB Version: 9.9.0.1495850 (R2020b) Update 1 was used to generate Figures 1, 2a, 2b, and 3.
2. The data and MATLAB codes required to generate Figures 2a, 2b, and 3.
3. A modified version of GROMACS-2021 that implements the solvers presented in the paper as an alternative to the default ones, SHAKE and P-LINCS. 
4. The required files to run the [Lysozyme in water simulation](https://www.mdtutorials.com/gmsozyme/) described in the paper (sections 6.2.2 and 6.2.3).

# Abstract

We analyze the convergence of quasi-Newton methods in exact and finite precision arithmetic using three different techniques. We derive an upper bound for the stagnation level and we show that any sufficiently exact quasi-Newton method will converge quadratically until stagnation. In the absence of sufficient accuracy, we are likely to retain rapid linear convergence. We confirm our analysis by computing square roots and solving bond constraint equations in the context of molecular dynamics. In particular, we apply both a symmetric variant and Forsgren's variant of the simplified Newton method. This work has implications for the implementation of quasi-Newton methods regardless of the scale of the calculation or the machine.

# Section 6.1: Computing square roots

The folder "simulations" contains a subfolder called "square-roots" with this content:

```
square-roots
├── fig
│   └── figure_1.pdf
└── matlab
    ├── figure_1.m
    └── newton_sqrt.m
```

The script figure_1 generates file figure_1.pdf using the MATLAB function newton_sqrt.
The script figure_1 assumes that current directory is square-roots/matlab.

# Section 6.2: Constrained molecular dynamics

We provide the required files to run the [Lysozyme in water simulation](https://manual.gromacs.org/documentation/2021/install-guide/index.html) as described in the paper. These files are located in the [simulations/lysozyme](simulations/lysozyme) folder.

We recommend using the [simulations/lysozyme/md.mdp](simulations/lysozyme/md.mdp) file to run the md simulation. In order to use our constraint solvers instead of the default ones, you must set the `ILVES_EXPERIMENT` environment variable to `1` `2` depending the experiment you want to reproduce. Set the variable to `1` to reproduce the experiments described in section 6.2.2 and set it to `2` to reproduce the experiments described in section 6.2.3. Also, to set the tolerance of the solver, change the `shake-tol` parameter of [md.mdp](simulations/lysozyme/md.mdp).

## GROMACS

The majority of the code added to GROMACS can be found in [Ilves.h](GROMACS/src/gromacs/mdlib/Ilves.h) and [Ilves.cpp](GROMACS/src/gromacs/mdlib/Ilves.cpp)

To install GROMACS, you can follow the [installation guide provided by GROMACS](https://manual.gromacs.org/documentation/2021/install-guide/index.html). **It is essential to generate a double-precision installation** in order to be able to use the constraint solvers with tiny tolerances. This is achieved by appending the flag `-DGMX_DOUBLE=on` to the CMake command.

### Installation example:
```
cd GROMACS
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_DOUBLE=on
make
sudo make install
source /usr/local/gromacs/bin/GMXRC
```

## Section 6.2.2: A quasi-Newton method with a fixed symmetric approximation of the Jacobian

The folder "simulations" contains a subfolder called "lysozyme" which includes this content

```
├── data
│   ├── normEk.txt
│   ├── reshis.txt
│   └── rk.txt
├── fig
│   ├── figure_2a.pdf
│   ├── figure_2b.pdf
│   └── figure_3.pdf
├── matlab
│   └── figure_2_and_3.m
```

The script figure_2_and_3 generates Figure 2a, Figure 2b and Figure 3 from the data stored in the folder data.
The data was generated with software included in this repository.

### How to run
```
cd simulations/lysozyme
# Set the desired tolerance by changing the `shake-tol` parameter of the `.mdp` file.
export ILVES_EXPERIMENT=1 # Execute the experiments described in section 6.2.2. Use our constraint solvers instead of the default ones.
export export NSCONSTRAINTS_OUTPUT=50 # Print every NSCONSTRAINTS_OUTPUT time-steps
# Generate the `.tpr` file
gmx_d grompp -f "md.mdp" -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -o "tpr.tpr"
# Run the simulation
gmx_d mdrun -ntmpi 1 -s "tpr.tpr" -noconfout -nsteps 50000 -g gromacs.log
```

This simulation will produce three files:
    - `rehis.txt`: The evolution of the maximum relative constraint violation. There is one line per time-step and one column per iteration.
    - `rk.txt`: The evolution of ||Ek||. There is one line per time-step and one column per iteration.
    - `normEk.txt`: The evolution of ||z - yk|| / ||z||. There is one line per time-step and one column per iteration.

These output files are used to generate Figure 2 and Figure 3.

## Section 6.2.3: Forsgren's variant of the simplified Newton method

### How to run accuracy simulation

```
cd simulations/lysozyme
# Set the desired tolerance by changing the `shake-tol` parameter of the `.mdp` file.
export ILVES_EXPERIMENT=2 # Execute the experiments described in section 6.2.3. Use our constraint solvers instead of the default ones.
export ILVES_PRINT=1 # Print accuracy results
export export NSCONSTRAINTS_OUTPUT=2 # Print every NSCONSTRAINTS_OUTPUT time-steps
# Generate the `.tpr` file
gmx_d grompp -f "md.mdp" -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -o "tpr.tpr"
# Run the simulation
gmx_d mdrun -ntmpi 1 -s "tpr.tpr" -noconfout -nsteps 50000 -g gromacs.log
```

The accuracy simulations will produce three files: `newton.txt`, `simpl.txt`, and `forsgren.txt`. These files will contain the accuracy results of Newton's, Simplified-Newton's, and Forsgren's methods, respectively.

The files contain seven columns and one line per time-step and constraint-algorithm iteration. Description of the columns:

- *time-step*: The time-step.
- *iteration*: The current iteration of the constraint-algorithm ($k$).
- *||lagr-newton||*: The 2-norm of the final approximation of the Lagrange multipliers obtained using Newton's method. 
- *tau*: The maximum relative constraint violation at iteration $k$.  
- *||lagr-correction-iter||*: The 2-norm of the correction of Lagrange multipliers at iteration $k$.
- *||current-lagr||*: The 2-norm of the approximation of the Lagrange multipliers at iteration $k$.
- *||lagr-newton - current_lagr||*: Normwise absolute error of the approximation of the Lagrange multipliers at iteration $k$ against the final approximation of the Lagrange multipliers obtained using Newton's method.
- *||lagr-newton - current-lagr||/||lagr-newton||*: Normwise relative error of the approximation of the Lagrange multipliers at iteration $k$ against the final approximation of the Lagrange multipliers obtained using Newton's method.

### How to run performance simulation

```
cd simulations/lysozyme
# Set the desired tolerance by changing the `shake-tol` parameter of the `.mdp` file.
export USE_ILVES=2 # Execute the experiments described in section 6.2.3. Use our constraint solvers instead of the default ones.
export ILVES_PRINT=0 # Do not print accuracy results
export export NSCONSTRAINTS_OUTPUT=1 # Measure time on each time-step
# Generate the `.tpr` file
gmx_d grompp -f "tpr.tpr" -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -o "tpr.tpr"
# Run the simulation
gmx_d mdrun -ntmpi 1 -s "tpr.tpr" -noconfout -nsteps 50000 -g gromacs.log
```

The performance simulations will produce one file named `times.txt`. The file contains one column per solver. The first line displays the total microseconds taken by each solver, and the second line the total number of iterations executed by each solver.

# How to cite us
TODO: When the manuscript has been accepted and the infomation is available.
