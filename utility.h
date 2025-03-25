#ifndef UTILITY_H
#define UTILITY_H

#define N 64
#define K 16
#define MAX 64
#define COUNT 10

void printMatrix(double *A);
void printVector(double *v);
void identity(double *A);
double norm(double *v);
void ata(double *A, double *ATA);
void productAb(double *A, double *b, double *c);
void productAB(double *A, double *B, double *C);
void qrMethod(double *A, double *V, double *W);
void qrDecomposition(double *A, double *Q, double *R);
void powerMethodN(double *A, double *V, double *W);
void powerMethod(double *A, double *v, double *w);
void powerMethodDeflation(double *A, double *v, double *w, double *B);
void leftSingularVectors(double *A, double *V, double *W, double *U);
void evdTopK(double *A, double *V, double *Ak);
void svdTopK(double *U, double *V, double *W, double *Ak);
void readFile(const char *filename, double *input);
void writeFile(const char *filename, double *output);

#endif
