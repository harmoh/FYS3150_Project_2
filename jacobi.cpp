#include <cmath>
#include <armadillo>
#include "jacobi.h"

using namespace arma;

// Main function for performing the Jacobi algorithm
void jacobi_method(int n)
{
    mat A = zeros(n - 1, n - 1);
    mat R(n - 1, n - 1);

    double rho_max, rho_min;
    vec rho(n);
    vec V(n - 1);
    vec d(n - 1);
    double e;
    rho_min = 0;
    rho_max = 5.0;

    double h = (rho_max - rho_min) / n;
    e = -1/(h*h);

    for(int i = 0; i < n - 1; i++)
    {
        rho(i) = rho(0) + i * h;
    }
    for(int i = 0; i < n - 1; i++)
    {
        V(i) = rho(i) * rho(i); // Non-interaction case
        d(i) = 2/(h*h) + V(i);

    }

    // Initialize A
    A(0,0) = d(0);
    for(int i = 1; i < n - 1; i++)
    {
        A(i,i) = d(i);
        A(i,i-1) = e;
        A(i-1,i) = e;
    }
    //cout << "A after:\n" << A << endl;
    cout << "rho:\n" << rho << endl;

    int k = 0;
    int l = 0;

    // Initializing eigenvalue matrix
    R.eye();

    double max_diagonal = max_offdiagonal(n, A, &k, &l);
    cout << "Max: " << max_diagonal << endl;

    double epsilon = 1e-8;
    double max_iterations = n * n * n;
    int iterations = 0;

    while(max_diagonal > epsilon && iterations < max_iterations)
    {
        max_diagonal = max_offdiagonal(n, A, &k, &l);

        rotate(n, A, R, k, l);

        iterations++;
    }

    //cout << "R:\n" << R << endl;
    cout << "Number of iterations: " << iterations << endl;

    vec eigenvalues(n);
    for(int i = 0; i < n - 1; i++)
    {
        eigenvalues(i) = R(i,i);
    }
    eigenvalues = sort(eigenvalues);

    cout << "eigenvalues:\n" << eigenvalues << endl;

    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, A);  // find eigenvalues/eigenvectors
    eigval = sort(eigval);
    cout << "eigval:\n" << eigval << endl;
    //cout << "eigvec:\n" << eigvec << endl;
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

void rotate(int n, mat &A, mat&R, int k, int l)
{
    double s, c;

    if(A(k,l) != 0.0)
    {
        double t, tau;
        tau = (A(k,k) - A(l,l)) / (2 * A(k,l));
        if(tau > 0)
        {
            t = 1.0 / (tau + sqrt(1 + tau * tau));
        }
        else
        {
            t = -1.0 / (-tau + sqrt(1 + tau * tau));
        }

        c = 1.0 / (sqrt(1 + t * t));
        s = c * t;
    }
    else
    {
        c = 1.0;
        s = 0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);

    // Changing the matrix elements with indices k and l
    A(k,l) = 0;
    A(l,k) = 0;
    A(k,k) = a_kk * c * c - 2 * A(k,l) * c * s + a_ll * s * s;
    A(l,l) = a_ll * c * c + 2 * A(k,l) * c * s + a_kk * s * s;

    for(int i = 0; i < n; i++)
    {
        if(i != k && i != l)
        {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik * c - A(i,l) * s;
            A(k,i) = A(i,k);
            A(i,l) = a_il * c + A(i,k) * s;
            A(l,i) = A(i,l);
        }
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = r_ik * c - r_il * s;
        R(i,l) = r_il * c + r_ik * s;
    }
}
