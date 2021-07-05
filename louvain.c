#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if !defined(ARRAY_SIZE)
    #define ARRAY_SIZE(x) (sizeof((x)) / sizeof((x)[0]))
#endif

float *array(int n) {
    return (float*) calloc(n, sizeof(float));
}

float **array2D(int n_rows, int n_cols) {
    int k = 0;
    float **matrix = malloc(n_rows * sizeof (float *));

    matrix[0] = malloc(n_cols * n_rows * sizeof (float) );
    for (k=1; k < n_rows; k++) {
        matrix[k] = matrix[0] + n_cols*k;
    }
    return matrix;
}

float *range(int n) {
    float *list = array(n);
    for(size_t i=0;i<n;i++) {
        list[i] = i;
    }
    return list;
}

float *arange(float min_val, float max_val, float delta) {
    int n = (int)ceil((max_val-min_val)/delta);
    printf("%d\n",n);
    float *list = array(n);
    for(int i=0;i<n;i++) {
        list[i] = min_val + i*delta;
    }
    return list;
}

float sum(float (**A)) {
    // Number of rows
    int n_rows = sizeof(A)/sizeof(A[0]);
    // Number of columns
    int n_cols  = sizeof(A)/sizeof(A[0][0]);
    // Sum elements
    float sum = 0;
    for(int i=0; i<n_rows; i++) {
        for(int j=0; j<n_cols; j++) {
            sum = sum + A[i][j];
        }
    }
    return sum;
}

float modularity_louvain_und(int n_nodes, float (**A), float gamma, int seed) {

    // Defining random seed
    srand(seed);

    // Number of nodes
    /* int n_nodes = sizeof(A)/sizeof(A[0]); */
    printf("n_nodes=%d\n", n_nodes);
    // Weight of edges
    float w_e   = sum(A);
    // Hierarchy index
    int h       = 0;
    // Hierarchical module assignments
    float *ci    = range(n_nodes);

    for(int i=0; i<n_nodes; i++) {
        printf("%f, ", ci[i]);
    }
    printf("\n");

    return 0.0;
}

int main() {

    float *a = arange(0,5,0.1);
    printf("%ld\n", ARRAY_SIZE(a));

    /* for(size_t i=0; i<50; i++){ */
    /*     printf("%f, ", a[i]); */
    /* } */
    /* printf("\n"); */


    float **A = array2D(5,5);
    for(int i=0;i<5;i++) {
        for(int j=0;j<5;j++){
            if(i==j) {
                A[i][j]=1;
            }
            else {
                A[i][j]=(i+1)*(j+1);
            }
        }
    }

    for(int i=0;i<5;i++) {
        for(int j=0;j<5;j++){
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }


    modularity_louvain_und(5, A, 1.0, 0);


    return 1;
}
