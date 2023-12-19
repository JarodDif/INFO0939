#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct _data{
    size_t N;
    double *vals;
}data_t;

#pragma omp declare mapper(data_t data) \
              map(data) \
              map(data.vals[0:data.N])

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int world_size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(world_size != 2){
        if(rank == 0){printf("Need exactly two processes\n")}
        MPI_Finalize();
        return 1;
    }

    data_t test;
    test.N = 10;
    if((test.vals = malloc(sizeof(double) * test.N)) == NULL){
        MPI_Abort(MPI_COMM_WORLD,MPI_ERR_IO);
    }

    for(int i=0; i < test.N; i++){
        test.vals[i] = (rank == 0)?2*i:0;
    }

    #pragma omp target enter data map(tofrom:test){
        double *ptr_to_vals = test.vals;
        #pragma omp target datause_device_ptr(ptr_to_vals){
            if(rank == 0){
                MPI_Send(ptr_to_vals, test.N, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            }else{
                MPI_Recv(ptr_to_vals, test.N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, null);
            }
        }
    }

    if(rank == 1){
        for(int i=0; i < test.N; i++){
            printf("%6.3lf ", test.vals[i]);
        }
        printf("\n");
    }

    free(test.vals);

    MPI_Finalize();

    return 0;
}