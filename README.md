# INFO0939 Project : Wave equation scheme (2023-2024)

## Students

* BERTRAND Tom (s201900)
* DIFFERDANGE Jarod (s202356)

## File structure

Each version of the program is in its own folder
* `default` contains the sequential version with the two asked modification
* `mpi` contains the pure MPI implementation
* `omp` contains the pure OMP implementation (without GPU acceleration)
* `mpi_omp` contains the MPI+OMP hybrid implementation
* `gpu` contains the GPU accelerated version using OMP pragmas

The given inputs are located in `example_inputs`. 
We also added the SLURM scripts as well as some other examples created by ourselves.

For standard compilation and running of the four first versions of the code, a Makefile has been created.
E.g. `make run_default` will run the sequential code on the given `param_3d.txt`

You will find the source code of the different tools given for the project (`data2gmsh`).

### Different versions

As you can see `mpi` and `mpi_omp` have a second version with the `_subdomain_save` suffix.
These versions of the code, instead of sending all their data to one rank so that it saves the outputs correctly, simply makes each process save their own data inside different output files, prefixed with their coordinates in the MPI Cartesian topology. To recombine the files, we coded `merge_output.py` (which uses numpy) to combine them in one file (only the first output specified in a parameter file is merged). Note that these versions of the code do not accept any other output type than `ALL`.

## Q1 : implement trilinear interpolation

In each version of the code, we added the function `trilinear_interpolation` and its usage in `interpolate_inputmaps`.
We coded the trilinear interpolation to get the $c$ and $\rho$ values for the simulation instead of the `closest_index` approach.

## Q2 : Experimental study of stability

The folder `default` contains a folder called `stability_test`. Inside we wrote a python script that runs the sequential code on a fixed $\delta$ and fixed $c$ while changing the value of $\Delta_t$.

We find a value close to $\frac{1}{\sqrt{3}}$ as the limit for the Courant number.

## Q3 : CPU bound or memory bound

The update of one pressure value takes around $10$ floating point operations and $10$ memory accesses.

The update of the velocity values at one index takes around $10$ floating point operations and $8$ memory accesses.

Considering that floating point operations take $6$ CPU cycles and each memory access takes $100-1000$ CPU cycles, we clearly are memory bound.

## Q4 : CPU cache bottlenecks

In each version of the code, the triple-loops have been reordered from `m -> n -> p` to `p -> n -> m` as the memory access used in the `INDEX3D` macro put values with consecutive values of `m` near each other. This drastically decreases the amount of cache misses.

## Q5 : MPI version

The folder `mpi` contains our MPI implementation of the code.

The values of $P_x, P_y, P_z$ are computed by `MPI_Dims_Create` and will create good values if the amount of processes are well-chosen.
With this we mean that if one were to chose a prime number (e.g. $5$), the code will separate the domain into $5$ domains only along one axis ($P_x = 5, P_y = P_z = 1$). 
But given a number with multiple prime factors (e.g. $8$), the code will separate the domain into $8$ domains, but always $2$ along each axis ($P_x = P_y = P_z = 2$).

This supposes that the domain is not to asymmetric, which corresponds to most of our cases. One could exhaustively find the best $P_x, P_y, P_z$ (which minimizes $T/p * (P_x/N_x + P_y/N_y + P_z/N_z)$) by exhaustive search.

## Q6 : OpenMP version

The folder `omp` contains our OpenMP implementation of the code.

After a small amount of tests using different values for the `collapse` clause, the value of $2$ performed marginally better than $3$ (the maximum) or no `collapse` clause at all.

## Q7 : Hybrid MPI + OpenMP version

The folder `mpi_omp` contains our hybrid implementation.

## Q8 : Performance measures and scalability analysis

Multiple slurm jobs have been created to generate reports for the different types of analysis on different cases.

For scalability in the folder `mpi` you have scripts and data relevant to the strong and weak scalability.

## Q9 : Single GPU acceleration 

The folder `gpu` contains our GPU accelerated version

This version performs significantly better than the OpenMP CPU version, even with $16$ processes.

## Q10 : Custom situations

We have created a speed input map that simulates the R7 classroom. We use this input map to simulate a teacher emitting a high frequency sound from the front of the classroom.

Similarly, we have also created a simulation of Young's double slit experiment and we clearly see the creation of Huygens' wavelets

## Optional enhancements

### 1 chili : theoretical study of the scheme

We have an analytical proof, using the Von Neumann stability criteria that Co (the Courant number) should be between $-\frac{1}{\sqrt{3}}$ and $\frac{1}{\sqrt{3}}$.

### 2 chilis : more sophisticated boundary conditions

We do not intend to do this point.

### 3 chilis : Multi-GPU version

We have two code versions that are prefixed by `mpi_gpu`.
1. One naive version that does not use `use_device_ptr`.
1. A better version that uses `use_device_ptr` and therefore is faster for inter-GPU communication.
