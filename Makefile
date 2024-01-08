# This Makefile is used to compile the codes and run them on the simple example

CC = gcc
MPICC = mpicc
CCLANG = clang
SCOREP = scorep
CFLAGS = -O3 -DNDEBUG
LDFLAGS = -lm
LDFLAGS_OMP = -lm -fopenmp
LDFLAGS_GPU = -lm -fopenmp --offload-arch=sm_80
SRC_DEFAULT = default/fdtd.c
HDR_DEFAULT = default/fdtd.h
OUT_DEFAULT = default/fdtd
SRC_MPI = mpi/fdtd.c
HDR_MPI = mpi/fdtd.h
OUT_MPI = mpi/fdtd
OUT_MPI_SCOREP = mpi/fdtd_instr
SRC_OMP = omp/fdtd.c
HDR_OMP = omp/fdtd.h
OUT_OMP = omp/fdtd
OUT_OMP_SCOREP = omp/fdtd_instr
SRC_HYBRID = mpi_omp/fdtd.c
HDR_HYBRID = mpi_omp/fdtd.h
OUT_HYBRID = mpi_omp/fdtd
OUT_HYBRID_SCOREP = mpi_omp/fdtd_instr
SRC_GPU = gpu/fdtd.c
HDR_GPU = gpu/fdtd.h
OUT_GPU = gpu/fdtd
SRC_MGPU = mpi_gpu_single_point/fdtd.c
HDR_MGPU = mpi_gpu_single_point/fdtd.h
OUT_MGPU = mpi_gpu_single_point/fdtd
EXAMPLE_FOLDER = example_inputs/simple3d/
PARAM_FILE = param_3d.txt

.PHONY: all run_default run_mpi run_omp run_hybrid clean

default: $(OUT_DEFAULT)
mpi: $(OUT_MPI)
mpi_scorep: $(OUT_MPI_SCOREP)
omp: $(OUT_OMP)
omp_scorep: $(OUT_OMP_SCOREP)
hybrid: $(OUT_HYBRID)
hybrid_scorep: $(OUT_HYBRID_SCOREP)
gpu: $(OUT_GPU)
mgpu: $(OUT_MGPU)

$(OUT_DEFAULT): $(SRC_DEFAULT) $(HDR_DEFAULT)
	$(CC) $(CFLAGS) -o $(OUT_DEFAULT) $(SRC_DEFAULT) $(LDFLAGS)

$(OUT_MPI): $(SRC_MPI) $(HDR_MPI)
	$(MPICC) $(CFLAGS) -o $(OUT_MPI) $(SRC_MPI) $(LDFLAGS)

$(OUT_MPI_SCOREP): $(SRC_MPI) $(HDR_MPI)
	$(SCOREP) $(MPICC) $(CFLAGS) -o $(OUT_MPI) $(SRC_MPI) $(LDFLAGS)

$(OUT_OMP) : $(SRC_OMP) $(HDR_OMP)
	$(CC) $(CFLAGS) -o $(OUT_OMP) $(SRC_OMP) $(LDFLAGS_OMP)

$(OUT_OMP_SCOREP) : $(SRC_OMP) $(HDR_OMP)
	$(SCOREP) $(CC) $(CFLAGS) -o $(OUT_OMP) $(SRC_OMP) $(LDFLAGS_OMP)

$(OUT_HYBRID): $(SRC_HYBRID) $(HDR_HYBRID)
	$(MPICC) $(CFLAGS) -o $(OUT_HYBRID) $(SRC_HYBRID) $(LDFLAGS_OMP)

$(OUT_HYBRID_SCOREP): $(SRC_HYBRID) $(HDR_HYBRID)
	$(SCOREP) $(MPICC) $(CFLAGS) -o $(OUT_HYBRID) $(SRC_HYBRID) $(LDFLAGS_OMP)

$(OUT_GPU) : $(SRC_GPU) $(HDR_GPU)
	$(CCLANG) $(CFLAGS) -o $(OUT_GPU) $(SRC_GPU) $(LDFLAGS_GPU)

$(OUT_MGPU) : $(SRC_MGPU) $(HDR_MGPU)
	$(MPICC) $(CFLAGS) -o $(OUT_MGPU) $(SRC_MGPU) $(LDFLAGS_GPU)

#runs the default version
run_default: $(OUT_DEFAULT)
	(cd $(EXAMPLE_FOLDER) && ../../$(OUT_DEFAULT) $(PARAM_FILE))

# runs the mpi version with NPROCESS processes
run_mpi: $(OUT_MPI)
ifndef NPROCESS
	$(error NPROCESS is not defined. Please provide it using 'make run_mpi NPROCESS=4')
endif
	(cd $(EXAMPLE_FOLDER) && mpirun -np $(NPROCESS) ../../$(OUT_MPI) $(PARAM_FILE))

# runs the omp version and asks for OMP_NUM_THREADS to be set
run_omp: $(OUT_OMP)
ifndef OMP_NUM_THREADS
	$(error OMP_NUM_THREADS is not defined. Please provide it using `make run_omp OMP_NUM_THREADS=4`)
endif
	(cd $(EXAMPLE_FOLDER) && ../../$(OUT_DEFAULT) $(PARAM_FILE))

# runs the hybrid version and asks for NPROCESS and OMP_NUM_THREADS to be set
run_hybrid: $(OUT_HYBRID)
ifndef NPROCESS
	$(error NPROCESS is not defined. Please provide it using 'make run_hybrid NPROCESS=4 OMP_NUM_THREADS=2')
endif
ifndef OMP_NUM_THREADS
	$(error OMP_NUM_THREADS is not defined. Please provide it using 'make run_hybrid NPROCESS=4 OMP_NUM_THREADS=2')
endif
	(cd $(EXAMPLE_FOLDER) && mpirun -np $(NPROCESS) ../../$(OUT_HYBRID) $(PARAM_FILE))


# removes the build programs
clean:
	rm -f $(OUT_DEFAULT) $(OUT_MPI) $(OUT_OMP) $(OUT_HYBRID) $(OUT_GPU) $(OUT_MGPU)
