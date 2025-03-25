#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utility.h"

int main(int argc, char *argv[])
{
    // Allocate memory
    double *input = (double *)malloc(COUNT * N * N * sizeof(double));
    double *outputEvdPm = (double *)malloc(COUNT * N * N * sizeof(double));
    double *outputSvdPm = (double *)malloc(COUNT * N * N * sizeof(double));
    double *outputEvdQrm = (double *)malloc(COUNT * N * N * sizeof(double));
    double *outputSvdQrm = (double *)malloc(COUNT * N * N * sizeof(double));
    double *ATA = (double *)malloc(N * N * sizeof(double));
    double *Ak = (double *)malloc(N * N * sizeof(double));
    double *U = (double *)malloc(N * N * sizeof(double));
    double *V = (double *)malloc(N * N * sizeof(double));
    double *W = (double *)malloc(N * N * sizeof(double));

    // Read input from disk
    const char *inputFilename = "input.bin";
    readFile(inputFilename, input);

    for (int i = 0; i < COUNT; i++)
    {
        // Find current image
        double *A = input + i * N * N;

        // Top k PC by EVD and PM
        ata(A, ATA);
        powerMethodN(ATA, V, W);
        evdTopK(A, V, Ak);
        memcpy(outputEvdPm + i * N * N, Ak, N * N * sizeof(double));

        // Top k PC by SVD and PM
        ata(A, ATA);
        powerMethodN(ATA, V, W);
        leftSingularVectors(A, V, W, U);
        svdTopK(U, V, W, Ak);
        memcpy(outputSvdPm + i * N * N, Ak, N * N * sizeof(double));

        // Top k PC by EVD and QRM
        ata(A, ATA);
        qrMethod(ATA, V, W);
        evdTopK(A, V, Ak);
        memcpy(outputEvdQrm + i * N * N, Ak, N * N * sizeof(double));

        // Top k PC by SVD and QRM
        ata(A, ATA);
        qrMethod(ATA, V, W);
        leftSingularVectors(A, V, W, U);
        svdTopK(U, V, W, Ak);
        memcpy(outputSvdQrm + i * N * N, Ak, N * N * sizeof(double));
    }

    // Write output to disk
    const char *outputEvdPmFilename = "outputEvdPm.bin";
    const char *outputSvdPmFilename = "outputSvdPm.bin";
    const char *outputEvdQrmFilename = "outputEvdQrm.bin";
    const char *outputSvdQrmFilename = "outputSvdQrm.bin";
    writeFile(outputEvdPmFilename, outputEvdPm);
    writeFile(outputSvdPmFilename, outputSvdPm);
    writeFile(outputEvdQrmFilename, outputEvdQrm);
    writeFile(outputSvdQrmFilename, outputSvdQrm);

    // Release memory
    free(outputEvdPm), free(outputSvdPm);
    free(outputEvdQrm), free(outputSvdQrm);
    free(input), free(ATA), free(Ak);
    free(U), free(V), free(W);
    return EXIT_SUCCESS;
}
