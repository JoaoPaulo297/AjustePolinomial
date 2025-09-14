#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define EPSILON 1e-18


// Troca duas linhas da matriz
void swapRows(double** A, int row1, int row2, int lenA) {
    for (int i = 0; i < 2 * lenA; i++) {  // Troca todas as colunas (incluindo a matriz identidade)
        double temp = A[row1][i];
        A[row1][i] = A[row2][i];
        A[row2][i] = temp;
    }
}

// Eliminação progressiva
void forwardElimination(double** A, int lenA) {
    for (int i = 0; i < lenA; i++) {
        // Encontrar a linha com o maior pivô absoluto na coluna i
        int maxRow = i;
        for (int k = i + 1; k < lenA; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Trocar linhas se necessário
        if (maxRow != i) {
            swapRows(A, i, maxRow, lenA);
        }

        // Verificar se o pivô é zero (matriz singular)
        if (fabs(A[i][i]) < EPSILON) {
            printf("Erro: matriz singular ou mal condicionada.\n");
            return;
        }

        // Eliminação Gaussiana
        for (int j = i + 1; j < lenA; j++) {
            double fator = A[j][i] / A[i][i];
            for (int k = i; k < 2 * lenA; k++) { // Atualiza todas as colunas
                A[j][k] -= A[i][k] * fator;
            }
        }
    }
}

// Eliminação reversa
void backwardElimination(double** A, int lenA) {
    for (int i = lenA - 1; i >= 0; i--) {
        // Verificar se o pivô é zero (matriz singular)
        if (fabs(A[i][i]) < EPSILON) {
            printf("Erro: matriz singular ou mal condicionada.\n");
            return;
        }

        // Normalizar a linha pelo pivô
        for (int j = 2 * lenA - 1; j >= i; j--) {
            A[i][j] /= A[i][i];
        }

        // Eliminação reversa
        for (int j = i - 1; j >= 0; j--) {
            double fator = A[j][i];
            for (int k = i; k < 2 * lenA; k++) {
                A[j][k] -= A[i][k] * fator;
            }
        }
    }
}

// Função para calcular a matriz inversa
double** inv(double** M, int dimension) {
    // Alocar memória para a matriz aumentada [M | I]
    double** A = malloc(dimension * sizeof(double*));
    for (int i = 0; i < dimension; i++) {
        A[i] = malloc(2 * dimension * sizeof(double));
    }

    // Inicializar a matriz aumentada [M | I]
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            A[i][j] = M[i][j]; // Copiar a matriz original
            A[i][j + dimension] = (i == j) ? 1.0 : 0.0; // Matriz identidade
        }
    }

    // Aplicar eliminação progressiva e reversa
    forwardElimination(A, dimension);
    backwardElimination(A, dimension);

    // Extrair a matriz inversa da parte direita da matriz aumentada
    double** inverseMatrix = malloc(dimension * sizeof(double*));
    for (int i = 0; i < dimension; i++) {
        inverseMatrix[i] = malloc(dimension * sizeof(double));
        for (int j = 0; j < dimension; j++) {
            inverseMatrix[i][j] = A[i][j + dimension];
        }
    }

    // Liberar memória da matriz aumentada
    for (int i = 0; i < dimension; i++) {
        free(A[i]);
    }
    free(A);

    return inverseMatrix;
}

void leituraDasCoordenadas(FILE *inputFile, double **Xs, double **Ys) {
    int lenArray;
    fscanf(inputFile, "%d", &lenArray);
    // Aloca memória para os arrays Xs e Ys
    *Xs = malloc(lenArray * sizeof(double));
    *Ys = malloc(lenArray * sizeof(double));

    // Verifica se a alocação de memória foi bem-sucedida
    if (*Xs == NULL || *Ys == NULL) {
        fprintf(stderr, "Erro ao alocar memória.\n");
        exit(1); // Encerra o programa em caso de erro
    }

    // Lê as coordenadas do arquivo
    for (int i = 0; i < lenArray; i++) {
        fscanf(inputFile, "%lf,%lf", &(*Xs)[i], &(*Ys)[i]);
    }
}

double** poly(double* X, int k, int lenX){
    double** matrix = malloc(lenX*sizeof(double*));
    for(int i = 0; i < lenX; i++){
        matrix[i] = malloc((k+1) * sizeof(double));
    }
    
    for(int i = 0; i < lenX; i++){
        for(int j = 0; j < k+1; j++){
            if(j != k)
               matrix[i][j] = pow(X[i], k-j);
            else
               matrix[i][j] = 1;
        }
    }
    
    return matrix;
    
}

void multiplicaMatrizes(double **A, double **B, double **C, int rowsA, int colsA, int colsB) {
    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
            C[i][j] = 0;
            for (int k = 0; k < colsA; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

double** transposta(double** A, int m, int n) {
    double** T = malloc(m * sizeof(double*));
    for(int i = 0; i < m; i++){
        T[i] = malloc(n * sizeof(double));
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            T[j][i] = A[i][j];  // Troca A[i][j] → T[j][i]
        }
    }
    return T;
}


void ajusta(double *X, double *Y, int n, int k, double *P, double *tempo, double *SSE) {
    clock_t inicio = clock();

    int cols = k + 1;

    // Gera a matriz de Vandermonde usando a função poly
    double **A = poly(X, k, n);

    // Alocação das matrizes
    double **At = malloc(cols * sizeof(double *));
    double **AtA = malloc(cols * sizeof(double *));
    double **AtA_inv = malloc(cols * sizeof(double *));
    double *AtY = malloc(cols * sizeof(double));

    // Verifica se a alocação de memória foi bem-sucedida
    if (A == NULL || At == NULL || AtA == NULL || AtA_inv == NULL || AtY == NULL) {
        fprintf(stderr, "Erro na alocação de memória!\n");
        exit(1);
    }

    // Aloca memória para as linhas das matrizes
    for (int i = 0; i < cols; i++) {
        At[i] = malloc(n * sizeof(double));
        AtA[i] = malloc(cols * sizeof(double));
        AtA_inv[i] = malloc(cols * sizeof(double));
        if (At[i] == NULL || AtA[i] == NULL || AtA_inv[i] == NULL) {
            fprintf(stderr, "Erro na alocação de memória!\n");
            exit(1);
        }
    }

    // Calcula A^T (transposta de A)
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < n; j++) {
            At[i][j] = A[j][i];
        }
    }

    // Calcula A^T * A
    multiplicaMatrizes(At, A, AtA, cols, n, cols);

    AtA_inv = inv(AtA, cols);

    // Calcula A^T * Y
    for (int i = 0; i < cols; i++) {
        AtY[i] = 0;
        for (int j = 0; j < n; j++) {
            AtY[i] += At[i][j] * Y[j];
        }
    }

    // Calcula P = (A^T * A)^(-1) * A^T * Y
    for (int i = 0; i < cols; i++) {
        P[i] = 0;
        for (int j = 0; j < cols; j++) {
            P[i] += AtA_inv[i][j] * AtY[j];
        }
    }

    // Calcula SSE
    double *Y_pred = malloc(n * sizeof(double));
    if (Y_pred == NULL) {
        fprintf(stderr, "Erro na alocação de memória!\n");
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        Y_pred[i] = 0;
        for (int j = 0; j < cols; j++) {
            Y_pred[i] += P[j] * A[i][j];
        }
    }

    *SSE = 0;
    for (int i = 0; i < n; i++) {
        *SSE += (Y_pred[i] - Y[i]) * (Y_pred[i] - Y[i]);
    }

    *tempo = ((double)(clock() - inicio)) / CLOCKS_PER_SEC;

    // Liberação de memória
    free(Y_pred);
    for (int i = 0; i < cols; i++) {
        free(At[i]);
        free(AtA[i]);
        free(AtA_inv[i]);
    }
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(At);
    free(AtA);
    free(AtA_inv);
    free(AtY);
}


int main(int argc, char *argv[]) {
    
    

    return 0;
}