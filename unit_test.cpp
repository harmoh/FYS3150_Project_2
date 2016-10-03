#include <armadillo>
#include "unit_test.h"
#include "jacobi.h"

using namespace std;
using namespace arma;

// Check if "max_offdiagonal()" finds the highest offdiagonal values
void test()
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
        cout << "Success! The max offdiagonal function works." << endl;
    }
    else
    {
        cout << "Not correct value." << endl;
    }
}

// Check that two double are equal
bool is_equal(double a, double b)
{
    double epsilon = 1e-8;
    return fabs(a - b) < epsilon;
}
