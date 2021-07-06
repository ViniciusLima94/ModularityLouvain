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

void copy(float *a_cp, float *a, int n) {
    for(int i=0; i<n; i++) {
        a_cp[i] = a[i];
    }
}

void copy2D(float **a_cp, float **a, int n_rows, int n_cols) {
    for(int i=0; i<n_rows; i++) {
        for(int j=0; j<n_cols; j++) {
            a_cp[i][j] = a[i][j];
        }
    }
}

float *append(float *a, int n, float val) {
    float *a_new = array(n+1);
    // Copy the old values to the appended array
    for(int i=0; i<n; i++) {
        a_new[i] = a[i];
    }
    // Append the new value
    a_new[n] = val;
    return a_new;
}

float max(float *a, int n) {
    float _max=0;
    for(int i=0; i<n; i++) {
        if(a[i]>_max) {
            _max=a[i];
        }
        else {
            continue;
        }
    }
    return _max;
}

int argmax(float *a, int n) {
    int _max    = 0;
    int _argmax = 0;
    for(int i=0; i<n; i++) {
        if(a[i]>_max) {
            _max    = a[i];
            _argmax = i;
        }
        else {
            continue;
        }
    }
    return _argmax;
}

void swap(float *xp, float *yp) {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// A function to implement bubble sort
// https://www.geeksforgeeks.org/bubble-sort/
float *sort(float *a, int n) {
    float *arr = array(n);
    copy(arr,a,n);
    int i, j;
    for (i = 0; i < n-1; i++) {
    // Last i elements are already in place
        for (j = 0; j < n-i-1; j++) {
            if (arr[j] > arr[j+1]) { swap(&arr[j], &arr[j+1]); }
        }
    }
    return arr;
}

int count_unique(float *a, int n) {
    int i,j,n_unq = 1;
    for(i=1;i<n;i++) {
        for(j=0;j<i;j++) {
            if(fabs(a[i]-a[j])<1e-10) {break;}
        }
        if(i==j) {n_unq++;}
    }
    return n_unq;
}

float *unique(float *a, int n, int sorted) {
    // Number of unique elements
    int N_UNQ = (int)count_unique(a,n);
    float *list_unq = array(N_UNQ);
    /* float list_unq[N_UNQ];  */
    list_unq[0] = a[0]; 

    int i,j,n_unq=1;
    for(i=1;i<n;i++) {
        for(j=0;j<i;j++) {
            if(fabs(a[i]-a[j])<1e-10) {break;}
        }
        if(i==j) {
            list_unq[n_unq]=a[i];
            n_unq++;
        }
    }
    if(sorted==0) {
        return list_unq;
    }
    else {
        return sort(list_unq,n_unq); 
    }
}

float *unique_reverse(float *a, int n) {
    // Get number of unique elements
    int N_UNQ = (int)count_unique(a,n);
    // Get unique elements not-sorted
    float *list_unq = unique(a,n,1);
    // Store the labels for each element
    float *rev_unq  = array(n);

    for(int i=0;i<n;i++) {
        for(int j=0;j<N_UNQ;j++) {
            if(a[i]==list_unq[j]) {rev_unq[i]=j;continue;}
        }
    }
    return rev_unq;
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
    float *list = array(n);
    for(int i=0;i<n;i++) {
        list[i] = min_val + i*delta;
    }
    return list;
}

float *get_row(float (**a), int row, int n_cols) {
    float *row_vals = array(n_cols);
    for(int i=0; i<n_cols; i++) {
        row_vals[i] = a[row][i];
    }
    return row_vals;
}

float *add_to_array(float *a, int n, float val) {
    float *a_sub = array(n);
    for(int i=0; i<n; i++) {
        a_sub[i] = a[i] + val;
    }
    return a_sub;
}

float *const_to_array(float *a, int n, float val) {
    float *a_prod = array(n);
    for(int i=0; i<n; i++) {
        a_prod[i] = a[i]*val;
    }
    return a_prod;
}

float **const_to_array2D(float (**a), int n, float val) {
    float **a_prod = array2D(n,n);
    for(int i=0; i<n; i++) {
        for(int j=0;j<n;j++) {
            a_prod[i][j] = a[i][j]*val;
        }
    }
    return a_prod;
}

float *add_array(float *a, float *b, int n) {
    float *c = array(n);
    for(int i=0; i<n; i++) {
        c[i] = a[i] + b[i];
    }
    return c;
}

float trace(float (**a), int n_rows) {
    float _trace = 0;
    for(int i=0; i<n_rows; i++) {
        _trace += a[i][i];
    }
    return _trace;
}

float *where(float *a, int n, float val) {
    int count = 0;
    float *indexes = array(1);
    for(int i=0; i<n; i++) {
        if(a[i]==val) {
            if(count==0) { indexes[0]=i; count++;}
            else { float *aux=append(indexes,n,i); indexes=aux;count++; }
        }
    }
    return indexes;
}

float *permutation(float *a, int n) {
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

float **multiply(float (**a) , float (**b), int N) {
    int i, j, k;
    float **res = array2D(N,N);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            res[i][j] = 0;
            for (k = 0; k < N; k++)
                res[i][j] += a[i][k] * b[k][j];
        }
    }
    return res;
}

float modularity_louvain_und(int n_nodes, float (**A), float gamma, int seed) {
    // Defining random seed
    srand(seed);
    // Weight of edges
    float w_e   = sum(A,n_nodes,n_nodes);
    // Hierarchy index
    int h       = 0;
    // Hierarchical module assignments
    float *ci   = range(n_nodes);
    // Hierarchical modularity values
    float q       = -1;
    // Nodes index
    float *n_i  = range(n_nodes);
    // 
    int n0      = n_nodes;

    while(1) {
        if(h>300) { printf("Entered an Infite Loop (E) - Aborted"); break; }
        // Nodes' degree
        float *k    = sum_alog_axis(A,n_nodes,n_nodes,0);
        // Copying values
        float *Km   = array(n_nodes); copy(Km,k,n_nodes);
        float **Knm = array2D(n_nodes, n_nodes); copy2D(Knm,A,n_nodes,n_nodes);
        // Initial module assignement
        float *m    = range(n_nodes);

        // Flag for within hierarchy search
        int flag = 1;
        int it   = 0;
        while(flag) {
            it++;
            if(it>1000) { printf("Entered an Infite Loop (F) - Aborted!"); break; }
            flag = 0;
            // Loop over nodes in random order
            float *n_ip = permutation(n_i, n_nodes);
            for(int i=0; i<n_nodes;i++) {
                int ii      = (int)n_ip[i];
                int ma      = (int)m[ii];
                // Algorithm condition
                float *dQ_1 = add_to_array(get_row(Knm,ii,n_nodes), n_nodes, -Knm[ii][ma]+A[ii][ii]);
                float *dQ_2 = const_to_array(add_to_array(Km,n_nodes,-Km[ma]+k[ii]),n_nodes,-gamma*k[ii]/w_e);
                float *dQ   = add_array(dQ_1,dQ_2,n_nodes);
                dQ[ma]      = 0;

                // Finding maximal modularity increae
                float max_dq = max(dQ,n_nodes);
                /* printf("%f ", max_dq); */
                // If maximal increase is positive
                if(max_dq>1e-10) {
                    // Take only one value
                    int j = argmax(dQ,n_nodes);
                    // Change modules degrees
                    for(int r=0;r<n_nodes;r++) {
                        Knm[r][j]  += A[r][ii];
                        Knm[r][ma] -= A[r][ii];
                    }
                    // Change module degrees
                    Km[j]  += k[ii];
                    Km[ma] -= k[ii];

                    // Reassign module
                    m[ii] = j;//j+1;
                    flag  = 1;
                }
            }
        }
        // New module assignments
        m = unique_reverse(m,n_nodes);
        /* m = add_to_array(m,n_nodes,1); */
        h++;
        float *ci_old = array(n_nodes); copy(ci_old,ci,n_nodes);
        ci            = array(n0);
        for(int i=0;i<n_nodes;i++) {
            for(int j=0;j<n_nodes;j++) {
                if(ci_old[j]==i) {
                    ci[j]=m[i];
                }
            }
        }
        // New number of modules
        n_nodes  = max(m,n_nodes);
        // New weighted matrix
        float **A1 = array2D(n_nodes,n_nodes);
        for(int i=0;i<n_nodes;i++) {
            for(int j=i;j<n_nodes;j++) {
                float wp = 0;
                // Pool weights of nodes in same module
                for(int r=0;r<n_nodes;r++) { if(m[r]==i && m[r]==j){ wp+=A[i][j]; } }
                A1[i][j] = wp;
                A1[j][i] = wp;
                printf("%f ", wp);
            }
        }
        copy2D(A,A1,n_nodes,n_nodes);

        // Compute modularity
        float q_old   = q;
        float **aux = const_to_array2D(A,n_nodes,1/w_e);
        q = trace(A,n_nodes)/w_e - gamma*sum(multiply(aux,aux,n_nodes),n_nodes,n_nodes);
        /* printf("%f ", q); */
        free(aux);
        if(q-q_old<1e10) {
            break;
        }
    }
    return q;
}

float readmatrix(int rows, int cols, float (**a), const char* filename) {

    FILE *pf;
    pf = fopen (filename, "-r");
    if (pf == NULL)
        return 0;

    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++)
            fscanf(pf, "%f", a[i] + j);
    }


    fclose (pf);
    return 1;
}


int main() {

    /* srand(1020); */

    int n     = 242;
    float **A = array2D(n,n);

    readmatrix(n,n,A,"matrix.txt");
    /* modularity_louvain_und(n, A, 1.0, 0); */
    /* float *k    = sum_alog_axis(A,n,n,0); */
    /* for(int r=0;r<n;r++) { printf("%f,", k[r]); } */
    /* printf("\n"); */

    /* for(int i=0;i<n;i++) { */
    /*     for(int j=0;j<n;j++) { */
    /*         printf("%d ", (int)A[i][j]); */
    /*     } */
    /* } */

    /* printf("%f ", trace(A,n)); */
    printf("%f ", modularity_louvain_und(n, A, 1.0, 10));

    //float *a = array(7);
    /* float a[8] = {1,1,2,5,5,3,3,6}; */
    /* int n      = ARRAY_SIZE(a); */
    /* float *indexes = where(a,n,5); */
    /* for(int i=0; i<2; i++) { */
    /*     printf("%d ", (int)indexes[i]); */
    /* } */
    /* printf("%d \n", n); */
    /* float *unq = unique_reverse(a,n); */
    /* printf("%d ", (int)count_unique(a,n)); */
    /* float *unq = unique(a,n,1); */
    /* float *unq = unique_reverse(a, n); */
    /* printf("\n"); */
    /* for(int i=0; i<n; i++) { */
    /*     printf("%d ", (int)unq[i]); */
    /* } */
    /* float *as = sort(a,n); */
    /* for(int i=0; i<8; i++) { */
    /*     printf("%d ", (int)as[i]); */
    /* } */


    /* int N_UNQ = (int)count_unique(a,n); */
    /* float list_unq[N_UNQ]; */
    /* list_unq[0] = a[0]; */
    /* int i,j,n_unq=1; */
    /* for(i=1;i<n;i++) { */
    /*     for(j=0;j<i;j++) { */
    /*         if(fabs(a[i]-a[j])<1e-10) {break;} */
    /*     } */
    /*     if(i==j) { */
    /*         [> float *list_unq = append(list_unq, n_unq, a[i]); <] */
    /*         list_unq[n_unq] = a[i]; */
    /*         printf("%f \n", list_unq[n_unq]); */
    /*         n_unq++; */
    /*     } */
    /* } */

    /* printf("\n"); */
    /* for(int i=0; i<n_unq; i++) { */
    /*     printf("%d ", (int)list_unq[i]); */
    /* } */


    /* count_unique(a,n); */
    /* float *a_unq = unique(a,n); */
    /* for(int i=0;i<n;i++){ */
    /*     printf("%d ", (int)a_unq[i]); */
    /* } */

    /* for(int i=0;i<n;i++){ */
    /*     printf("%d ", (int)a[i]); */
    /* } */
    /* printf("\n"); */
    /* for(int i=0;i<n;i++){ */
    /*     printf("%d ", (int)unq[i]); */
    /* } */

    /* printf("%d \n", count_unique(a,ARRAY_SIZE(a))); */

    /* a[5]=10; */
    /* float _max = max(a,10); */
    /* printf("%f \n", _max); */
    /* int _argmax = argmax(a,10); */
    /* printf("%d \n", _argmax); */
    /* float *a   = arange(0,5,1); */

    /* float *pa  = permutation(a,5); */

    /* float *aa  = append(a,5,11); */
    
    /* float *acp = array(6); */
    /* copy(acp,aa,6); */

    /* float *as = add_to_array(a,5,-3); */

    /* float *ap = const_to_array(a,5,-3); */

    /* float *c  = add_array(a,ap,5); */

    /* for(size_t i=0; i<5; i++){ */
    /*     printf("%f, ", a[i]); */
    /* } */
    /* printf("\n"); */

    /* for(size_t i=0; i<5; i++){ */
    /*     printf("%f, ", pa[i]); */
    /* } */
    /* printf("\n"); */

    /* for(size_t i=0; i<6; i++){ */
    /*     printf("%f, ", aa[i]); */
    /* } */
    /* printf("\n"); */

    /* for(size_t i=0; i<6; i++){ */
    /*     printf("%f, ", acp[i]); */
    /* } */
    /* printf("\n"); */

    /* for(size_t i=0; i<5; i++){ */
    /*     printf("%f, ", as[i]); */
    /* } */
    /* printf("\n"); */

    /* for(size_t i=0; i<5; i++){ */
    /*     printf("%f, ", ap[i]); */
    /* } */
    /* printf("\n"); */

    /* for(size_t i=0; i<5; i++){ */
    /*     printf("%f, ", c[i]); */
    /* } */
    /* printf("\n"); */
    /* printf("\n"); */

    /* for(int i=0;i<5;i++) { */
    /*     for(int j=0;j<5;j++){ */
    /*         if(i==j) { */
    /*             A[i][j]=1; */
    /*         } */
    /*         else { */
    /*             A[i][j]=(i+1)*(j+1); */
    /*         } */
    /*     } */
    /* } */

    /* float **Acp = array2D(5,5); */
    /* copy2D(Acp,A,5,5); */

    /* for(int i=0;i<5;i++) { */
    /*     for(int j=0;j<5;j++){ */
    /*         printf("%f ", A[i][j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */
    /* printf("\n"); */

    /* for(int i=0;i<5;i++) { */
    /*     for(int j=0;j<5;j++){ */
    /*         printf("%f ", Acp[i][j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */
    /* printf("\n"); */

    /* float *row_vals = get_row(A,3,5); */

    /* for(size_t i=0; i<5; i++){ */
    /*     printf("%f, ", row_vals[i]); */
    /* } */
    /* printf("\n"); */

    /* float *sum = sum_alog_axis(A, 5, 5, 1); */

    /* for(int i=0;i<5;i++) { */
    /*     printf("%f ",sum[i]); */
    /* } */
    /*     printf("\n"); */

    /* modularity_louvain_und(5, A, 1.0, 0); */

    return 1;
}
