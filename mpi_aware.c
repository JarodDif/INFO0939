#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct _data{
    size_t N;
    double *vals;
    double **neighbors;
}data_t;

#pragma omp declare mapper(data_t data) \
              map(data) \
              map(data.vals[0:data.N]) \
              map(data.neighbors[0:1]) \
              map(data.neighbors[0][0:data.N])

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
    test.N = 10;
    if((test.vals = malloc(sizeof(double) * test.N)) == NULL &&
        (test.neighbors = malloc(sizeof(double*) * 6)) == NULL  &&
        (test.neighbors[0] = malloc(sizeof(double) * test.N)) == NULL){
        MPI_Abort(MPI_COMM_WORLD,MPI_ERR_IO);
    }

    for(int i=0; i < test.N; i++){
        test.vals[i] = (rank == 0)?2*i:0;
        test.neighbors[0][i] = (rank == 0)?2*i:0;
    }

    #pragma omp target data map(tofrom:test)
    {
        double *ptr_to_vals = test.vals;
        double *ptr_to_n = test.neighbors[0];
        #pragma omp target data use_device_ptr(ptr_to_vals, ptr_to_n)
        {
            if(rank == 0){
                MPI_Send(ptr_to_vals, test.N, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
                MPI_Send(ptr_to_n, test.N, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
            }else{
                MPI_Recv(ptr_to_vals, test.N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, NULL);
                MPI_Recv(ptr_to_n, test.N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, NULL);
            }
        }
    }

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