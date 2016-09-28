#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;
ofstream ofile_d;

// Main function for finding eigenvalues using Armadillo for the interacting case
void interacting()
{
    // Initial variables
    double rho_min = 0;
    double rho_max = 5.0;
    vec w_r(4);
    w_r(0) = 0.01;
    w_r(1) = 0.5;
    w_r(2) = 1.0;
    w_r(3) = 5.0;

    //    clock_t start_arma, finish_arma;

    int n = 100;
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

    // Loop over different w_r values
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < n-1; j++)
        {
            V(j) = w_r(i)*w_r(i) * rho(j+1)*rho(j+1) + 1/rho(j+1); // Interacting case
            d(j) = 2/(h*h) + V(j);
        }

        // Initialize A
        A(n-2,n-2) = d(n-2);
        for(int i = 0; i < n-2; i++)
        {
            A(i,i) = d(i);
            A(i,i+1) = e;
            A(i+1,i) = e;
        }

        //    start_arma = clock();
        vec eigenval;
        mat eigenvec;

        eig_sym(eigenval, eigenvec, A);  // find eigenvalues/eigenvectors
        eigenval = sort(eigenval);

        // Normalizing eigenvectors with a simple trapezoidal rule:
        double norm1 = 0;
        double norm2 = 0;
        double norm3 = 0;
        vec V1 = eigenvec.col(0);
        vec V2 = eigenvec.col(1);
        vec V3 = eigenvec.col(2);

        vec U_square_normed1(n-1);
        vec U_square_normed2(n-1);
        vec U_square_normed3(n-1);

        for (int i = 0; i < n-1; i++) {
            double V1_i2 = V1(i)*V1(i);
            double V2_i2 = V2(i)*V2(i);
            double V3_i2 = V3(i)*V3(i);

            norm1 += V1_i2;
            norm2 += V2_i2;
            norm3 += V3_i2;

            U_square_normed1(i) = V1_i2;
            U_square_normed2(i) = V2_i2;
            U_square_normed3(i) = V3_i2;
        }
        norm1 *= h;
        norm2 *= h;
        norm3 *= h;

        U_square_normed1 /= norm1;
        U_square_normed2 /= norm2;
        U_square_normed3 /= norm3;

        // Writing all the normalized results to a txt file and plot them
        // with a python script.
        string outfilename = "Results_interacting_case_";
        stringstream w_r_string;
        w_r_string << w_r(i);
        outfilename.append(w_r_string.str());
        outfilename.append(".txt");
        cout << "w_r(" << i << "): " << w_r(i) << endl;
        cout << "outfilename: " << outfilename << endl;
        ofile_d.open(outfilename);
        ofile_d << setiosflags(ios::showpoint | ios::uppercase);
        ofile_d << " rho:               u1:                 u2:                 u3:" << endl;
        // Loop over all n producing table of time used:
        for (int i=0; i<n-1; i++) {
            double U_val1 = U_square_normed1(i);
            double U_val2 = U_square_normed2(i);
            double U_val3 = U_square_normed3(i);
            ofile_d << setw(10) << setprecision(8) << rho(i+1);
            ofile_d << setw(20) << setprecision(8) << U_val1;
            ofile_d << setw(20) << setprecision(8) << U_val2;
            ofile_d << setw(20) << setprecision(8) << U_val3 << endl;
        }
        ofile_d.close();
    }
}
