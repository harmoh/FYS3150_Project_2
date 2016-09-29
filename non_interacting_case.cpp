#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include "non_interacting_case.h"

using namespace std;
using namespace arma;
ofstream ofile_plot;

// Main function for calculating eigenvalues and eigenvectors for the non-interacting case
void non_interacting()
{
    // Initial variables
    double rho_min = 0;
    double rho_max = 5.0;

    int n = 350;
    mat A = zeros(n-1, n-1);
    mat R(n-1, n-1);

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

    // Initializing eigenvalue matrix
    R.eye();

    clock_t start_arma, finish_arma;

    start_arma = clock();
    vec eigenval;
    mat eigenvec;

    eig_sym(eigenval, eigenvec, A);  // find eigenvalues/eigenvectors
    eigenval = sort(eigenval);
    finish_arma = clock();
    double time_arma = (finish_arma - start_arma)/(double)CLOCKS_PER_SEC;

    cout << "Time: " << time_arma << endl;

    // Normalizing eigenvectors:
    double norm1 = 0;
    double norm2 = 0;
    double norm3 = 0;
    vec eigenvec1 = eigenvec.col(0);
    vec eigenvec2 = eigenvec.col(1);
    vec eigenvec3 = eigenvec.col(2);

    vec u_square_norm1(n-1);
    vec u_square_norm2(n-1);
    vec u_square_norm3(n-1);

    for (int j = 0; j < n-1; j++) {
        double eigenvec1_squared = eigenvec1(j)*eigenvec1(j);
        double eigenvec2_squared = eigenvec2(j)*eigenvec2(j);
        double eigenvec3_squared = eigenvec3(j)*eigenvec3(j);

        norm1 += eigenvec1_squared;
        norm2 += eigenvec2_squared;
        norm3 += eigenvec3_squared;

        u_square_norm1(j) = eigenvec1_squared;
        u_square_norm2(j) = eigenvec2_squared;
        u_square_norm3(j) = eigenvec3_squared;
    }
    norm1 *= h;
    norm2 *= h;
    norm3 *= h;

    u_square_norm1 /= norm1;
    u_square_norm2 /= norm2;
    u_square_norm3 /= norm3;

    // Writing all the normalized results to a txt file. One file for each w_r value.
    string outfilename = "Results_non_interacting_case";
    //outfilename.append(to_string(n));
    outfilename.append(".txt");
    ofile_plot.open(outfilename);
    ofile_plot << setiosflags(ios::showpoint | ios::uppercase);
    ofile_plot << " rho:               u1:                 u2:                 u3:" << endl;
    // Loop over all n producing table of time used:
    for (int j = 0; j < n-1; j++) {
        double u_val1 = u_square_norm1(j);
        double u_val2 = u_square_norm2(j);
        double u_val3 = u_square_norm3(j);
        ofile_plot << setw(10) << setprecision(8) << rho(j+1);
        ofile_plot << setw(20) << setprecision(8) << u_val1;
        ofile_plot << setw(20) << setprecision(8) << u_val2;
        ofile_plot << setw(20) << setprecision(8) << u_val3 << endl;
    }
    ofile_plot.close();
}
