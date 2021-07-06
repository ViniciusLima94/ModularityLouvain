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

float *permutation(float *a, int n)
{
    
    // Allocate and initialize permutated array
    float *perm_array = array(n);
    for(int i=0; i<n; i++) {
        perm_array[i] = a[i];
    }
    if(n > 1) {
        for (int i=0; i<n-1; i++)
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t    = perm_array[j];
          perm_array[j] = perm_array[i];
          perm_array[i] = t;
        }
    }
    return perm_array;
}

float sum(float (**A), int n_rows, int n_cols) {
    // Sum elements
    float _sum = 0;
    for(int i=0; i<n_rows; i++) {
        for(int j=0; j<n_cols; j++) {
            _sum = _sum + A[i][j];
        }
    }
    return _sum;
}

float *sum_alog_axis(float (**A), int n_rows, int n_cols, int axis) {
    if(axis==0) {
        float *_sum = array(n_cols);
        for(int i=0; i<n_cols; i++) {
            for(int j=0; j<n_rows; j++) {
                _sum[i] = _sum[i] + A[j][i];
            }
        }
        return _sum;
    }
    else {
        float *_sum = array(n_rows);
        for(int i=0; i<n_rows; i++) {
            for(int j=0; j<n_cols; j++) {
                _sum[i] = _sum[i] + A[i][j];
            }
        }
        return _sum;
    }
}

float modularity_louvain_und(int n_nodes, float (**A), float gamma, int seed) {

    // Defining random seed
    srand(seed);

    // Number of nodes
    /* int n_nodes = sizeof(A)/sizeof(A[0]); */
    printf("n_nodes=%d\n", n_nodes);
    // Weight of edges
    float w_e   = sum(A, n_nodes, n_nodes);
    // Hierarchy index
    int h       = 0;
    // Hierarchical module assignments
    float *ci   = range(n_nodes);
    // Nodes index
    float *n_i  = range(n_nodes);

    /* for(int i=0; i<n_nodes; i++) { */
    /*     printf("%f, ", ci[i]); */
    /* } */
    /* printf("\n"); */

    while(1) {
        if(h>300) {
            printf("Entered an Infite Loop (E) - Aborted");
            break;
        }
        float *k_o = sum_alog_axis(A,n_nodes,n_nodes,1);
        float *k_i = sum_alog_axis(A,n_nodes,n_nodes,0);
        // Initial module assignement
        float *m = range(n_nodes);

        // Flag for within hierarchy search
        int flag = 1;
        int it   = 0;
        while(flag) {
            it++;
            if(it>1000) {
                printf("Entered an Infite Loop (F) - Aborted!");
                break;
            }
            flag = 0;
        }
        // Loop over nodes in random order
        float *n_ip = permutation(n_i, n_nodes);
        for(int i=0; i<n_nodes;i++) {
            float ma = m[(int)n_ip[i]];
            // Algorithm condition
            float dQ;
            for(int i=0; i<n_nodes; i++) {
                for(int j=0; j<n_nodes; j++) {
                    (A[i,j]-A[i,ma]+A[i,i])-gamma*k[i]*(k[i
                }
            }
        }

    }

    return 0.0;
}

int main() {

    float *a = arange(0,5,1);
    printf("%ld\n", ARRAY_SIZE(a));

    float *pa = permutation(a,5);

    for(size_t i=0; i<5; i++){
        printf("%f, ", a[i]);
    }
    printf("\n");

    for(size_t i=0; i<5; i++){
        printf("%f, ", pa[i]);
    }
    printf("\n");


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

    /* for(int i=0;i<5;i++) { */
    /*     for(int j=0;j<5;j++){ */
    /*         printf("%f ", A[i][j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */

    /* float *sum = sum_alog_axis(A, 5, 5, 1); */

    /* for(int i=0;i<5;i++) { */
    /*     printf("%f ",sum[i]); */
    /* } */
    /*     printf("\n"); */

    modularity_louvain_und(5, A, 1.0, 0);


    return 1;
}
