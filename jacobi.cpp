#include <cmath>
#include <armadillo>
#include "jacobi.h"

using namespace arma;

// Main function for performing the Jacobi algorithm
void jacobi_method(int n)
{
    mat A = zeros(n-1, n-1);
    mat R(n-1, n-1);

    double rho_max, rho_min;
    vec rho(n+1);
    vec V(n-1);
    vec d(n-1);
    double e;
    rho_min = 0;
    rho_max = 5.0;

    double h = (rho_max - rho_min) / n;
    e = -1.0/(h*h);

    for(int i = 0; i < n+1; i++)
    {
        rho(i) = rho(0) + i * h;
    }
    for(int i = 0; i < n-1; i++)
    {
        V(i) = rho(i + 1) * rho(i + 1); // Non-interaction case
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

    // Initializing eigenvalue matrix
    R.eye();

    int k = 0;
    int l = 0;

    double epsilon = 1.0e-8;
    double max_iterations = n * n * n;
    int iterations = 0;

    double max_diagonal = max_offdiagonal(n-1, A, &k, &l);
    while(max_diagonal > epsilon && iterations < max_iterations)
    {
        max_diagonal = max_offdiagonal(n-1, A, &k, &l);

        rotate(n-1, A, R, k, l);

        iterations++;
    }

    cout << "\nNumber of iterations: " << iterations << endl;

    vec eigenvalues(n-1);
    for(int i = 0; i < n-1; i++)
    {
        eigenvalues(i) = A(i,i);
    }
    eigenvalues = sort(eigenvalues);

    cout << "Eigenvalues:\t";
    for(int i = 0; i < 3; i++)
    {
        cout << eigenvalues(i) << "\t";
    }
    cout << endl;

    mat B = A;
    eigenvalues_arma(B);
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
    double s, c, t, tau;;

    if(A(k,l) != 0.0)
    {
        tau = (A(l,l) - A(k,k)) / (2 * A(k,l));
        if(tau >= 0)
        {
            t = 1.0 / (fabs(tau) + sqrt(1 + tau * tau));
        }
        else
        {
            t = -1.0 / (fabs(tau) + sqrt(1 + tau * tau));
        }

        c = 1.0 / (sqrt(1.0 + t * t));
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
    A(k,k) = a_kk * c * c - 2 * A(k,l) * c * s + a_ll * s * s;
    A(l,l) = a_ll * c * c + 2 * A(k,l) * c * s + a_kk * s * s;
    A(k,l) = 0;
    A(l,k) = 0;

    for(int i = 0; i < n; i++)
    {
        if(i != k && i != l)
        {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik * c - a_il * s;
            A(k,i) = A(i,k);
            A(i,l) = a_il * c + a_ik * s;
            A(l,i) = A(i,l);
        }
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = r_ik * c - r_il * s;
        R(i,l) = r_il * c + r_ik * s;
    }
}

void eigenvalues_arma(mat &B)
{
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, B);  // find eigenvalues/eigenvectors
    eigval = sort(eigval);
}
