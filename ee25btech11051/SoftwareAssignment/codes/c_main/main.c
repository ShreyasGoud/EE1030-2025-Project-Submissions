#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "stb_image.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"functions.h"


int main(){

    
    int columns, rows, channels;

    unsigned char *gray_img = stbi_load("globe.jpg", &columns, &rows, &channels, 1);
    if (gray_img == NULL) {
        printf("Error: Could not load image\n");
        return 1;
    }
    
    int**matrix = (int**)calloc(rows, sizeof(int*));
    for(int i=0; i<rows; i++){
        *(matrix+i) = (int*)calloc(columns, sizeof(int));
    }

    double**constmatrix = (double**)calloc(rows, sizeof(double*));
    for(int i=0; i<rows; i++){
        *(constmatrix+i) = (double*)calloc(columns, sizeof(double));
    }

    

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < columns; x++) {
            matrix[y][x] = (int)gray_img[y * columns + x];
        }
    }

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < columns; x++) {
            constmatrix[y][x] =(double) matrix[y][x];
        }
    }
    
    int n=columns;
    int m=rows;

    for(int i=0; i<5; i++){
        for(int j=0; j<5; j++){
            printf("%d", matrix[i][j]);
        }
    }
    printf("\n");

     size_t bytes = 
         (size_t)m*(size_t)n*sizeof(double) * 2 +
        (size_t)n*(size_t)n*sizeof(double) * 2 +
         (size_t)m*(size_t)n*sizeof(double);

     printf("Approx memory required: %.1f MB\n", bytes / 1048576.0);

     if (bytes > 80000000ULL) { 
         printf("ERROR: Image too large for available memory on this system.\n");
        exit(1);
    }

    double**V = (double**)calloc(n, sizeof(double*));
    for(int i=0; i<n; i++){
        V[i]=(double*)calloc(n, sizeof(double));
    }

    double**U = (double**)calloc(m, sizeof(double*));
    for(int i=0; i<m; i++){
        U[i]=(double*)calloc(n, sizeof(double));
    }
    
    double**S = (double**)calloc(n, sizeof(double*));
    for(int i=0; i<n; i++){
        S[i]=(double*)calloc(n, sizeof(double));
    }

    double **Awork = (double**)calloc(m, sizeof(double*));
    for(int i=0; i<m; i++){
        Awork[i] = (double*)calloc(n, sizeof(double));
        for(int j=0; j<n; j++)
            Awork[i][j] = (double)matrix[i][j];
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                *(*(V+i)+j)=1;
            }else{
                *(*(V+i)+j)=0;
            }
        }
    }

    double toler = 1e-6;
    int max_sweeps = 1500;
    int k = 3;

    jacobi(Awork, max_sweeps, toler, m, n, V, S, U);
    sort(n, U, V, S, m, n);

    double **US = (double**)calloc(m, sizeof(double*));
    for (int i = 0; i < m; i++)
        US[i] = (double*)calloc(k, sizeof(double));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            US[i][j] = U[i][j] * S[j][j];
        }
    }

    double **A_k = (double**)calloc(m, sizeof(double*));
    for (int i = 0; i < m; i++)
        A_k[i] = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int r = 0; r < k; r++) {
                sum += US[i][r] * V[j][r];
            }
            A_k[i][j] = sum;
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = A_k[i][j];
            if (v < 0) v = 0;
            if (v > 255) v = 255;
            matrix[i][j] = (int)(v + 0.5);
        }
    }

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < columns; x++) {
            constmatrix[y][x] = constmatrix[y][x]-matrix[y][x];
        }
    }

    double norm = frobeniusNorm(constmatrix, rows, columns);
    
    unsigned char *rec_img = (unsigned char*)calloc((size_t)m*(size_t)n, sizeof(unsigned char));
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            rec_img[i*n + j] = (unsigned char)matrix[i][j];
        }
    }

    
    char filename[64];

    snprintf(filename, sizeof(filename), "globe-%d.jpg", k);

    int quality = 95; 
    if(stbi_write_jpg(filename, n, m, 1, rec_img, quality)) {
        
        // You should also update your printf to be correct!
        printf("Saved reconstructed image as %s\n", filename);
    
    } else {
        printf("Error saving JPG image.\n");
    }

    printf("Frobenius norm = %lf", norm);

    printf("Cleaning up memory...\n");
    stbi_image_free(gray_img);
    free(rec_img);
    
    free_double_matrix(constmatrix, m);
    free_int_matrix(matrix, m);
    free_double_matrix(V, n);
    free_double_matrix(U, m);
    free_double_matrix(S, n);
    free_double_matrix(Awork, m);
    free_double_matrix(US, m);
    free_double_matrix(A_k, m);

    printf("\nReconstruction complete. Saved as %s\n", filename);


    return 0;
}