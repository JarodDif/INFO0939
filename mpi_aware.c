#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct _data{
    size_t N;
    double *vals;
    double **neighbors;
}data_t;

typedef struct _sim{
    data_t *data;
}sim_t;

#pragma omp declare mapper(data_t data) \
              map(data) \
              map(data.vals[0:data.N]) \
              map(data.neighbors[0:1]) \
              map(data.neighbors[0][0:data.N])

#pragma omp declare mapper(sim_t sim) \
            map(sim) \
            map(sim.data[0:1])

void doTransfer(sim_t *sim, int rank){

    double *ptr_to_vals = sim->data->vals;
    double *ptr_to_n = sim->data->neighbors[0];
    #pragma omp target data use_device_ptr(ptr_to_vals, ptr_to_n)
    {
        if(rank == 0){
            MPI_Send(ptr_to_vals, sim->data->N, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            MPI_Send(ptr_to_n, sim->data->N, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
        }else{
            MPI_Recv(ptr_to_vals, sim->data->N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, NULL);
            MPI_Recv(ptr_to_n, sim->data->N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, NULL);
        }
    }

}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int world_size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(world_size != 2){
        if(rank == 0){printf("Need exactly two processes\n");}
        MPI_Finalize();
        return 1;
    }

    data_t test;
    sim_t sim;
    sim.data = &test;
    test.N = 10;
    if((test.vals = malloc(sizeof(double) * test.N)) == NULL ||
        (test.neighbors = malloc(sizeof(double*) * 6)) == NULL  ||
        (test.neighbors[0] = malloc(sizeof(double) * test.N)) == NULL){
        MPI_Abort(MPI_COMM_WORLD,MPI_ERR_IO);
    }

    printf("%d : malloced\n", rank);

    for(int i=0; i < test.N; i++){
        test.vals[i] = (rank == 0)?2*i:0;
        test.neighbors[0][i] = (rank == 0)?2*i:0;
    }

    printf("%d : filled data \n", rank);

    #pragma omp target data map(tofrom:sim)
    {
        doTransfer(&sim, rank); 
    }  

    printf("%d : Finished transfer\n", rank);

    if(rank == 1){
        for(int i=0; i < test.N; i++){
            printf("%6.3lf ", test.vals[i]);
        }
        printf("\n");
        for(int i=0; i < test.N; i++){
            printf("%6.3lf ", test.neighbors[0][i]);
        }
        printf("\n");
    }

    free(test.vals);
    free(test.neighbors[0]); free(test.neighbors);

    MPI_Finalize();

    return 0;
}