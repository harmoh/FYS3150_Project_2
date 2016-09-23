#include <cmath>
#include <armadillo>
#include "jacobi.h"

using namespace arma;

// Main function for performing the Jacobi algorithm
void jacobi_method(int n)
{
    mat A(n,n);
    mat R(n,n);

    vec ro(n);
    vec V(n);
    vec d(n);
    double e;
    ro(0) = 0;
    ro(n - 1) = 5.0;

    double h = (ro(n - 1) - ro(0)) / n;
    e = -1/(h*h);

    for(int i = 1; i < n; i++)
    {
        ro(i) = ro(0) + i * h;
    }
    for(int i = 0; i < n; i++)
    {
        V(i) = ro(i) * ro(i); // Non-interaction case
        d(i) = 2/(h*h) + V(i);

    }

    // Initialize A
    A(0,0) = d(0);
    for(int i = 1; i < n; i++)
    {
        A(i,i) = d(i);
        A(i,i-1) = e;
        A(i-1,i) = e;
    }
    cout << "A after:\n" << A << endl;

    int k = 0;
    int l = 0;

    // Initializing eigenvalue matrix
    R.eye();

    double max_diagonal = max_offdiagonal(n, A, &k, &l);
    cout << "Max: " << max_diagonal << endl;

    double epsilon = 10e-8;
    double max_iterations = n * n * n;
    int iterations = 0;

    while(max_diagonal > epsilon && iterations < max_iterations)
    {
        max_diagonal = max_offdiagonal(n, A, &k, &l);


        iterations++;
    }

    cout << "Number of iterations: " << iterations << endl;
}

double max_offdiagonal(int n, mat A, int *k, int *l)
{
    double max = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = i + 1; j < n; j++)
        {
            if(fabs(A(i,j)) > max)
            {
                max = fabs(A(i,j));
                *k = i;
                *l = j;
            }
        }
    }
    return max;
}

void rotate()
{
    double s, c;

    if(A(k,l) != 0)
    {
        double t, tau;
        tau = (A(k,k) - A(l,l)) / (2 * A(k,l));
    }
}
