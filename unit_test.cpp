#include <armadillo>
#include "unit_test.h"
#include "jacobi.h"

using namespace std;
using namespace arma;

// Check if "max_offdiagonal()" finds the highest offdiagonal values
void test_max_off()
{
    // Unit test
    // Create and initialze a nxn tridiagonal matrix with random values
    int n = 6;
    int k = 0;
    int l = 0;
    mat A = zeros(n-1, n-1);
    double d_n = rand();
    A(n-2, n-2) = d_n;
    double max_temp = 0;
    for(int i = 0; i < n-2; i++)
    {
        double d_i = rand();
        double e_ij = rand();
        double e_ji = rand();
        A(i,i) = d_i;
        A(i,i+1) = e_ij;
        A(i+1,i) = e_ji;

        if(fabs(e_ij) > max_temp)
        {
            max_temp = fabs(e_ij);
        }
        if(fabs(e_ji) > max_temp)
        {
            max_temp = fabs(e_ji);
        }
    }
    double max = max_offdiagonal(n-1, A, &k, &l);

    cout << "Unit test 1: ";
    if(is_equal(max, max_temp))
    {
        cout << "Success! The max offdiagonal function returns the maximum value." << endl;
    }
    else
    {
        cout << "Not correct value from the offdiagonal function." << endl;
    }
}

// Initialize same matrix A as in the Jacobi method and check that diagonal values are not zero
void test_non_empty()
{
    // Initial variables
    double rho_min = 0;
    double rho_max = 5.0;

    int n = 6;
    mat A = zeros(n-1, n-1);

    vec rho(n+1);
    rho(0) = rho_min;
    vec V(n-1);
    vec d(n-1);

    double h = (rho_max - rho_min) / n;
    double e = -1.0/(h*h);

    for(int i = 0; i < n+1; i++)
    {
        rho(i) = rho(0) + i * h;
    }
    for(int i = 0; i < n-1; i++)
    {
        V(i) = rho(i+1) * rho(i+1); // Non-interaction case
        d(i) = 2/(h*h) + V(i);
    }

    // Initialize A
    A(n-2,n-2) = d(n-2);
    for(int i = 0; i < n-2; i++)
    {
        A(i,i) = d(i);
        A(i,i+1) = e;
        A(i+1,i) = e;
    }

    cout << "Unit test 2: ";
    bool non_empty = true;
    for(int i = 0; i < n-2; i++)
    {
        if(is_equal(A(i,i), 0))
        {
            non_empty = false;
        }
    }
    if(non_empty)
    {
        cout << "Success! The diagonal values of matrix A are not empty." << endl;
    }
    else
    {
        cout << "Error! The diagonal values of matrix A are not empty." << endl;
    }
}

// Check that two double values are equal
bool is_equal(double a, double b)
{
    double epsilon = 1e-8;
    return fabs(a - b) < epsilon;
}
