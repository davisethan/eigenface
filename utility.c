#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include "utility.h"

/**
 * @brief Print matrix
 * @param A Matrix to print
 */
void printMatrix(double *A)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f ", A[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * @brief Print vector
 * @param v Vector to print
 */
void printVector(double *v)
{
    for (int i = 0; i < N; i++)
    {
        printf("%f ", v[i]);
    }
    printf("\n");
}

/**
 * @brief Make identity matrix
 * @param A Matrix to make identity
 */
void identity(double *A)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i * N + j] = i == j ? 1.0 : 0.0;
        }
    }
}

/**
 * @brief Euclidean norm of vector
 * @param v Vector
 * @return Euclidean norm
 */
double norm(double *v)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

/**
 * @brief Matrix transpose matrix
 * @param A Matrix
 * @param ATA Matrix transpose matrix
 */
void ata(double *A, double *ATA)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < N; k++)
            {
                sum += A[k * N + i] * A[k * N + j];
            }
            ATA[i * N + j] = sum;
        }
    }
}

/**
 * @brief Matrix vector multiplication
 * @param A Matrix
 * @param b Vector
 * @param c Product vector
 */
void productAb(double *A, double *b, double *c)
{
    for (int i = 0; i < N; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < N; j++)
        {
            sum += A[i * N + j] * b[j];
        }
        c[i] = sum;
    }
}

/**
 * @brief Matrix matrix multiplication
 * @param A Left matrix
 * @param B Right matrix
 * @param C Product matrix
 */
void productAB(double *A, double *B, double *C)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < N; k++)
            {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

/**
 * @brief QR method
 * @param A Matrix
 * @param V Eigenvectors matrix
 * @param W Eigenvalues matrix
 */
void qrMethod(double *A, double *V, double *W)
{
    // Allocate memory
    double *Q = (double *)malloc(N * N * sizeof(double));
    double *R = (double *)malloc(N * N * sizeof(double));
    double *C = (double *)malloc(N * N * sizeof(double));

    // Initialize memory
    identity(V);

    for (int i = 0; i < MAX; i++)
    {
        qrDecomposition(A, Q, R);

        // Update matrix
        productAB(R, Q, C);
        memcpy(A, C, N * N * sizeof(double));

        // Update eigenvectors
        productAB(V, Q, C);
        memcpy(V, C, N * N * sizeof(double));
    }

    // Update eigenvalues
    memcpy(W, A, N * N * sizeof(double));

    // Release memory
    free(Q), free(R), free(C);
}

/**
 * @brief QR decomposition
 * @param A Matrix
 * @param Q Orthogonal matrix
 * @param R Upper triangular matrix
 */
void qrDecomposition(double *A, double *Q, double *R)
{
    // LAPACK orthogonal matrix Q
    memcpy(Q, A, N * N * sizeof(double));
    double u[N];
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, N, N, Q, N, u);
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, N, N, N, Q, N, u);

    // LAPACK upper triangular matrix R
    memcpy(R, A, N * N * sizeof(double));
    double v[N];
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, N, N, R, N, v);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            R[i * N + j] = 0.0;
        }
    }
}

/**
 * @brief Power method iterations
 * @param A Matrix
 * @param V Eigenvectors matrix
 * @param W Eigenvalues matrix
 */
void powerMethodN(double *A, double *V, double *W)
{
    // Allocate memory
    double *B = (double *)malloc(N * N * sizeof(double));
    double *v = (double *)malloc(N * sizeof(double));
    double w = 0.0;

    for (int i = 0; i < N; i++)
    {
        powerMethod(A, v, &w);

        // Update eigenvectors
        for (int j = 0; j < N; j++)
        {
            V[j * N + i] = v[j];
        }

        // Update eigenvalues
        W[i * N + i] = w;

        powerMethodDeflation(A, v, &w, B);
        memcpy(A, B, N * N * sizeof(double));
    }

    // Release memory
    free(B), free(v);
}

/**
 * @brief Power method
 * @param A Matrix
 * @param v Dominant eigenvector
 * @param w Dominant eigenvalue
 */
void powerMethod(double *A, double *v, double *w)
{
    // Allocate memory
    double *u = (double *)malloc(N * sizeof(double));

    // Initialize memory
    for (int i = 0; i < N; i++)
    {
        v[i] = 1.0;
    }

    for (int i = 0; i < MAX; i++)
    {
        // Update dominant eigenvector
        productAb(A, v, u);
        memcpy(v, u, N * sizeof(double));

        // Normalize dominant eigenvector
        double nv = norm(v);
        for (int j = 0; j < N; j++)
        {
            v[j] = v[j] / nv;
        }
    }

    // Find paired dominant eigenvalue
    productAb(A, v, u);
    double a = 0.0, b = 0.0;
    for (int i = 0; i < N; i++)
    {
        a += v[i] * u[i];
        b += v[i] * v[i];
    }
    *w = a / b;

    // Release memory
    free(u);
}

/**
 * @brief Power method deflation
 * @param A Matrix to deflate
 * @param v Dominant eigenvector
 * @param w Dominant eigenvalue
 * @param B Deflated matrix
 */
void powerMethodDeflation(double *A, double *v, double *w, double *B)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            B[i * N + j] = A[i * N + j] - *w * v[i] * v[j];
        }
    }
}

/**
 * @brief Left singular vectors
 * @param A Matrix
 * @param V Right singular vectors
 * @param W Singular values squared
 * @param U Left singular vectors
 */
void leftSingularVectors(double *A, double *V, double *W, double *U)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < N; k++)
            {
                sum += A[i * N + k] * V[k * N + j];
            }
            U[i * N + j] = sum / sqrt(W[j * N + j]);
        }
    }
}

/**
 * @brief Top k principal components of matrix
 * @param A Matrix
 * @param V Eigenvectors matrix
 * @param Ak Matrix with top k principal components
 */
void evdTopK(double *A, double *V, double *Ak)
{
    double *B = (double *)malloc(N * N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double sum = 0.0;
            for (int l = 0; l < K; l++)
            {
                sum += V[i * N + l] * V[j * N + l];
            }
            B[i * N + j] = sum;
        }
    }
    productAB(A, B, Ak);
    free(B);
}

/**
 * @brief Top k principal components of matrix
 * @param U Left singular vectors
 * @param V Right singular vectors
 * @param W Singular values squared
 * @param Ak Matrix with top k principal components
 */
void svdTopK(double *U, double *V, double *W, double *Ak)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double sum = 0.0;
            for (int l = 0; l < K; l++)
            {
                sum += U[i * N + l] * sqrt(W[l * N + l]) * V[j * N + l];
            }
            Ak[i * N + j] = sum;
        }
    }
}

/**
 * @brief Read images from disk
 * @param filename Images filename
 * @param input Memory space of images
 */
void readFile(const char *filename, double *input)
{
    FILE *file = fopen(filename, "rb");
    fread(input, sizeof(double), COUNT * N * N, file);
    fclose(file);
}

/**
 * @brief Write images to disk
 * @param filename Images filename
 * @param output Memory space of images
 */
void writeFile(const char *filename, double *output)
{
    FILE *file = fopen(filename, "wb");
    fwrite(output, sizeof(double), COUNT * N * N, file);
    fclose(file);
}
