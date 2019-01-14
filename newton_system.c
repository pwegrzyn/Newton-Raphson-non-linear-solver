#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define EPSILON 0.00001
#define N 3

double jacobi_matrix[N * N];
double X[N];
double X_prim[N];
double B[N];
double X_last[N];

double F1(double v[N]);
double F2(double v[N]);
double F3(double v[N]);
void gaussian_elimination(double *a, double *b, double *x, int n);
void newton_raphson_system(double[N]);
double derivative(double (*f)(double *), double vect[N], int j);
void matrix_vector_mul(double *matrix, int a, int b, double *vect, int n, double *res);
double euclidean_norm(double[N], double[N]);
void copy_vect(double *a, double *b, int n);
void generate_jacobi_matrix(double vect[N]);

double F1(double v[N]) {
    return pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2) - 1;
}

double F2(double v[N]) {
    return 2 * pow(v[0], 2) + pow(v[1], 2) - 4 * v[2];
}

double F3(double v[N]) {
    return 3 * pow(v[0], 2) - 4 * v[1] + pow(v[2], 2);
}

// Calculate the derivative of the function with respect to the j-th argument in point vect
double derivative(double (*f)(double *), double vect[N], int j) {
    double new_vect_1[N], new_vect_2[N];
    for(int i = 0; i < N; i++) {
        new_vect_1[i] = vect[i];
        new_vect_2[i] = vect[i];
    }
    const double delta = 1.0e-6;
    new_vect_1[j] = vect[j] - delta;
    new_vect_2[j] = vect[j] + delta;
    double y1 = f(new_vect_1);
    double y2 = f(new_vect_2);
    double res = (y2 - y1) / (2 * delta);
    return res;
}

// Solve the given non-linear system of equations using Newtons method and init initial solution
void newton_raphson_system(double init[N]) {
    double zeros[N];
    for(int i = 0; i < N; i++) {
        X[i] = init[i];
        zeros[i] = 0;
    }
    // if the initial vector is [0, 0, ..., 0] then change it a bit so that the jacobi matrix will become invertible 
    if(euclidean_norm(X, zeros) < EPSILON) {
        for(int i = 0; i < N; i++) {
            X[i] += (double)rand() / RAND_MAX;
        }
    }
    int counter = 0;
    do {
        copy_vect(X, X_last, N);
        generate_jacobi_matrix(X);
        matrix_vector_mul(jacobi_matrix, N, N, X, N, B);
        X_prim[0] = B[0] - F1(X);
        X_prim[1] = B[1] - F2(X);
        X_prim[2] = B[2] - F3(X);
        gaussian_elimination(jacobi_matrix, X_prim, X, N);
        counter++;
    } while (euclidean_norm(X, X_last) >= EPSILON);
    printf("Newton-Raphson iterations: %d, epsilon: %f\n", counter, EPSILON);
    printf("Found solution: [ ");
    for(int i = 0; i < N; i++) {
        printf("%f ", X[i]);
    }
    printf("]\n");
}

void generate_jacobi_matrix(double vect[N]) {
    for(int i = 0; i < N; i++) {
        jacobi_matrix[0 * N + i] = derivative(F1, vect, i);
        jacobi_matrix[1 * N + i] = derivative(F2, vect, i);
        jacobi_matrix[2 * N + i] = derivative(F3, vect, i);
    }
}

void copy_vect(double *a, double *b, int n) {
    for(int i = 0; i < n; i++) {
        b[i] = a[i];
    }
}

void gaussian_elimination(double *a, double *b, double *x, int n) {
    int column, row, diagonal, max_pivot_row, j;
    double max_pivot, tmp;
    for (diagonal = 0; diagonal < n; diagonal++)
    {
        max_pivot_row = diagonal;
        max_pivot = *(a + (diagonal * n + diagonal)); // i,ith element of the matrix
        for (row = diagonal + 1; row < n; row++)
        {
            tmp = fabs(*(a + (row * n + diagonal)));
            if (tmp > max_pivot)
            {
                max_pivot_row = row;
                max_pivot = tmp;
            }
        }

        if (diagonal != max_pivot_row)
        {
            for (int k = 0; k < n; k++)
            {
                double *tmp_pointer1 = a + (diagonal * n + k);
                double *tmp_pointer2 = a + (max_pivot_row * n + k);
                tmp = *tmp_pointer1;
                *tmp_pointer1 = *tmp_pointer2;
                *tmp_pointer2 = tmp;
            }
            tmp = b[diagonal];
            b[diagonal] = b[max_pivot_row];
            b[max_pivot_row] = tmp;
        }

        for (row = diagonal + 1; row < n; row++)
        {
            tmp = *(a + (row * n + diagonal)) / *(a + (diagonal * n + diagonal));
            for (column = diagonal + 1; column < n; column++)
            {
                *(a + (row * n + column)) -= tmp * *(a + (diagonal * n + column));
            }
            *(a + (row * n + diagonal)) = 0;
            b[row] -= tmp * b[diagonal];
        }
    }

    for (row = n - 1; row >= 0; row--)
    {
        tmp = b[row];
        for (j = n - 1; j > row; j--)
        {
            tmp -= x[j] * *(a + (row * n + j));
        }
        x[row] = tmp / *(a + (row * n + row));
    }
}

void matrix_vector_mul(double *matrix, int a, int b, double *vect, int n, double *res)
{
    if (b != n)
    {
        perror("Shapes do not match");
        exit(1);
    }
    for (int i = 0; i < a; i++)
    {
        res[i] = 0;
        for (int j = 0; j < b; j++)
        {
            res[i] += matrix[i * b + j] * vect[j];
        }
    }
}

double euclidean_norm(double x_1[N], double x_2[N]) {
    double sum = 0;
    for(int i = 0; i < N; i++) {
        sum += pow(x_1[i] - x_2[i], 2);
    }
    return sqrt(sum);
}

int main(void) 
{
    srand(time(NULL));
    double x_init[] = {-15000, 8, 1242.9};
    newton_raphson_system(x_init);

}