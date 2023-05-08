# [PAPER] Newton's method revisited: How accurate do we have to be?

This repository contains the code used to generate the data for Figures 1, 2, 3, 4, 5, and 6 of the paper titled [Newton's method revisited: How accurate do we have to be?](TODO). The repository contains the MATLAB codes necessary to compute square roots and generate Figure 1. The repository contains a modified version of GROMACS-2021 that implements the solvers presented in the paper as an alternative to the default ones, SHAKE and P-LINCS. It also contains the required files to run the [Lysozyme in water simulation](https://www.mdtutorials.com/gmx/lysozyme/) described in the paper.

# Abstract

We analyze the convergence of quasi-Newton methods in exact and finite precision arithmetic using three different techniques. We derive an upper bound for the stagnation level and we show that any sufficiently exact quasi-Newton method will converge quadratically until stagnation. In the absence of sufficient accuracy, we are likely to retain rapid linear convergence. We confirm our analysis by computing square roots and solving bond constraint equations in the context of molecular dynamics. In particular, we apply both a symmetric variant and Forsgren's variant of the simplified Newton method. This work has implications for the implementation of quasi-Newton methods regardless of the scale of the calculation or the machine.

# Computing square roots.

The folder "simulations" contains a subfolder called "square-roots" with this content:
square-roots
├── fig
│   └── figure_1.pdf
└── matlab
    ├── figure_1.m
    └── newton_sqrt.m

The script figure_1 generates file figure_1.pdf using the MATLAB function newton_sqrt.
The script figure_1 assumes that current directory is square-roots/matlab.

MATLAB is a commercial product and the property of MathWorks, hence it cannot be included in this repository.
MATLAB Version: 9.9.0.1495850 (R2020b) Update 1 was used to generate the figure.


# Modified GROMACS

The majority of the code added to GROMACS can be found in [Ilves.h](GROMACS/src/gromacs/mdlib/Ilves.h) and [Ilves.cpp](GROMACS/src/gromacs/mdlib/Ilves.cpp)

# Installation of GROMACS

To install GROMACS, you can follow the [installation guide provided by GROMACS](https://manual.gromacs.org/documentation/2021/install-guide/index.html). **It is essential to generate a double-precision installation** in order to be able to use the constraint solvers with tiny tolerances. This is achieved by appending the flag `-DGMX_DOUBLE=on` to the CMake command.

## Installation Example:
```
cd GROMACS
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_DOUBLE=on
make
sudo make install
source /usr/local/gromacs/bin/GMXRC
```

# Run Simulations

We provide the required files to run the [Lysozyme in water simulation](https://manual.gromacs.org/documentation/2021/install-guide/index.html) as described in the paper. These files are located in the [simulations/lysozyme](simulations/lysozyme) folder.

We recommend using the [simulations/lysozyme/md.mdp](simulations/lysozyme/md.mdp) file to run the md simulation. In order to use our constraint solvers instead of the default ones, you must set the `USE_ILVES` environment variable to `1`. Also, to set the tolerance of the solver, change the `shake-tol` parameter of [md.mdp](simulations/lysozyme/md.mdp).

## Run Accuracy Simulation

The accuracy simulations will produce three files: `newton.txt`, `quasi.txt`, and `forsgren.txt`. These files will contain the accuracy results of Newton's, Quasi-Newton's, and Forsgren's methods, respectively.

The files contain seven columns and one line per time-step and constraint-algorithm iteration. Description of the columns:

- *time-step*: The time-step.
- *iteration*: The current iteration of the constraint-algorithm ($k$).
- *||lagr-newton||*: The 2-norm of the final approximation of the Lagrange multipliers obtained using Newton's method. 
- *tau*: The maximum relative constraint violation at iteration $k$.  
- *||lagr-correction-iter||*: The 2-norm of the correction of Lagrange multipliers at iteration $k$.
- *||current-lagr||*: The 2-norm of the approximation of the Lagrange multipliers at iteration $k$.
- *||lagr-newton - current_lagr||*: Normwise absolute error of the approximation of the Lagrange multipliers at iteration $k$ against the final approximation of the Lagrange multipliers obtained using Newton's method.
- *||lagr-newton - current-lagr||/||lagr-newton||*: Normwise relative error of the approximation of the Lagrange multipliers at iteration $k$ against the final approximation of the Lagrange multipliers obtained using Newton's method.

### How to Run
```
cd simulations/lysozyme
# Set the desired tolerance by changing the `shake-tol` parameter of the `.mdp` file.
export USE_ILVES=1 # Use our constraint solvers instead of the default ones.
export ILVES_PRINT=1 # Print accuracy results
export export NSCONSTRAINTS_OUTPUT=2 # Print every NSCONSTRAINTS_OUTPUT time-steps
# Generate the `.tpr` file
gmx_d grompp -f "md.mdp" -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -o "tpr.tpr"
# Run the simulation
gmx_d mdrun -ntmpi 1 -s "tpr.tpr" -noconfout -nsteps 50000 -g gromacs.log
```

## Run Performance Simulation

The performance simulations will produce one file named `times.txt`. The file contains one column per solver. The first line displays the total microseconds taken by each solver, and the second line the total number of iterations executed by each solver.

### How to Run
```
cd simulations/lysozyme
# Set the desired tolerance by changing the `shake-tol` parameter of the `.mdp` file.
export USE_ILVES=1 # Use our constraint solvers instead of the default ones.
export ILVES_PRINT=0 # Do not print accuracy results
export export NSCONSTRAINTS_OUTPUT=1 # Measure time on each time-step
# Generate the `.tpr` file
gmx_d grompp -f "tpr.tpr" -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -o "tpr.tpr"
# Run the simulation
gmx_d mdrun -ntmpi 1 -s "tpr.tpr" -noconfout -nsteps 50000 -g gromacs.log
```

# How to Cite Us
TODO: When the manuscript has been accepted and the infomation is available.
