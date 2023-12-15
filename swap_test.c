#include <stdio.h>

#define N 10

void swap_pointers(int **a, int **b) {
    int *temp = *a;
    *a = *b;
    *b = temp;
}

int main() {
    int *host_a, *host_b;

    // Allocate memory on the host
    host_a = (int *)malloc(N * sizeof(int));
    host_b = (int *)malloc(N * sizeof(int));

    // Initialize data on the host
    for (int i = 0; i < N; i++) {
        host_a[i] = i;
        host_b[i] = 2 * i;
    }

    // Map pointers to device memory
    #pragma omp target enter data map(to: host_a[0:N], host_b[0:N])

    // Swap pointers on the device
    #pragma omp target teams distribute parallel for
    for (int i = 0; i < N; i++) {
        swap_pointers(&host_a, &host_b);
    }

    // Map pointers back to the host
    #pragma omp target exit data map(from: host_a[0:N], host_b[0:N])

    // Verify the result on the host
    printf("After swapping:\n");
    printf("host_a: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", host_a[i]);
    }
    printf("\n");

    printf("host_b: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", host_b[i]);
    }
    printf("\n");

    // Free allocated memory
    free(host_a);
    free(host_b);

    return 0;
}
