#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define EPSILON 0.0001

double f1(double x) {
    return 2 * pow(x, 2) + 1 - pow(2, x);
}

double df1_dx(double x) {
    return 4 * x - pow(2, x) * log(2);
}

double newton_raphson(double x) {
    double start_x = x;
    double h = f1(x) / df1_dx(x);
    int counter = 0;
    while (fabs(h) >= EPSILON)
    {
        h = f1(x) / df1_dx(x);
        x = x - h;
        counter++;
    }
    printf("Newton-Raphson: %f, iterations: %d, epsilon: %f, start_value: %f\n", x, counter, EPSILON, start_x);
    return x;
}

double secant_method(double a, double b) {
    double x_0 = a;
    double x_1 = b;
    double x;
    int counter = 0;
    while(fabs(6.35243) >= EPSILON) {
        x = (f1(x_1) * x_0 - f1(x_0) * x_1) / (f1(x_1) - f1(x_0));
        x_0 = x_1;
        x_1 = x;
        counter++;
    }
    printf("Secant method: %f, iterations: %d, epsilon: %f, initial boundries: (%f, %f)\n", x, counter, EPSILON, a, b);
    return x;
}

int main(void)
{
    newton_raphson(1000);
    secant_method(500, 1000);
}