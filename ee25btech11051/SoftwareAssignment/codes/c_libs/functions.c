#include<stdio.h>

double signum(double k){ //signum function
    if(k>0.0){
        return 1;
    }else if(k<0.0){
        return -1;
    }else return 0;
}

double frobeniusNorm(double **mat, int rows, int cols) {
    double sumSq = 0.0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            sumSq += mat[i][j] * mat[i][j]; // Sum of squares of elements
        }
    }
    return sqrt(sumSq); // Return the square root of the sum of squares
}

double** transpose(int n, int m, double ** matrix){ //finds the transpose of a matrix
    double**matrix2 = (double**)calloc(m, sizeof(double*));

    for(int i=0; i<m; i++){
        matrix2[i]=(double*)calloc(n, sizeof(double));
    }

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            matrix2[i][j]= matrix[j][i];
        }
    }
    return matrix2;
}

//The below jacobi's function gets the values for U, V and S given A.

void jacobi(double **A, int max_sweeps, double toler, int m, int n, double **V, double **S, double **U){
    for(int sweep = 0; sweep < max_sweeps; sweep++) {
        int converged = 1;
        for(int p = 0; p <= n - 2; p++) {
            for(int q = p + 1; q < n; q++) {

                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for(int i = 0; i < m; i++) {
                    double ap = A[i][p];
                    double aq = A[i][q];
                    alpha += ap * ap;
                    beta  += aq * aq;
                    gamma += ap * aq;
                }

                if (fabs(gamma) <= toler * sqrt(alpha * beta + 1e-30))
                    continue;

                converged = 0;
                double tau = (beta - alpha) / (2.0 * gamma);
                double t   = signum(tau) / (fabs(tau) + sqrt(1.0 + tau * tau));
                double c   = 1.0 / sqrt(1.0 + t * t);
                double s   = c * t;

                for(int i = 0; i < m; i++) {
                    double aip = A[i][p];
                    double aiq = A[i][q];
                    double newp = c * aip - s * aiq;
                    double newq = s * aip + c * aiq;
                    A[i][p] = newp;
                    A[i][q] = newq;
                }

                for(int i = 0; i < n; i++) {
                    double vip = V[i][p];
                    double viq = V[i][q];
                    double newp = c * vip - s * viq;
                    double newq = s * vip + c * viq;
                    V[i][p] = newp;
                    V[i][q] = newq;
                }
            }
        }
        if(converged)
            break;
    }

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            S[i][j] = 0.0;
    for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
            U[i][j] = 0.0;

    for(int j = 0; j < n; j++) {
        double sigma_sq = 0.0;
        for(int i = 0; i < m; i++)
            sigma_sq += A[i][j] * A[i][j];
        double sigma = sqrt(sigma_sq);
        S[j][j] = sigma;
        if(sigma > toler) {
            for(int i = 0; i < m; i++)
                U[i][j] = A[i][j] / sigma;
        }
    }
}

//A bubble sort that sorts the singular values and the eigenvalues of U and V.

void sort(int k, double**U, double**V, double**S, int m, int n){

    double *array = (double*)calloc(k, sizeof(double));
    for(int i=0; i<k; i++)
        array[i] = S[i][i];

    double *Uarray = (double*)calloc(m, sizeof(double));
    double *Varray = (double*)calloc(n, sizeof(double));

    for (int end = k-1; end >= 0; end--) {
        int swapper = 0;
        for (int j = 0; j < end; j++) {
            if (array[j] < array[j+1]) {
                swapper++;
                double temp = array[j];
                array[j] = array[j+1];
                array[j+1] = temp;

                for(int i=0; i<m; i++){
                    Uarray[i] = U[i][j];
                    U[i][j] = U[i][j+1];
                    U[i][j+1] = Uarray[i];
                }
                for(int i=0; i<n; i++){
                    Varray[i] = V[i][j];
                    V[i][j] = V[i][j+1];
                    V[i][j+1] = Varray[i];
                }
            }
        }
        if (swapper == 0) break;
    }

    for(int i=0; i<k; i++)
        S[i][i] = array[i];

    free(array);
    free(Uarray);
    free(Varray);
}

void free_double_matrix(double **matrix, int rows) {
    if (!matrix) return;
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void free_int_matrix(int **matrix, int rows) {
    if (!matrix) return;
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}