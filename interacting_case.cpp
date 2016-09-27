#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;
ofstream ofile_interacting;

// Main function for performing the Jacobi algorithm
void interacting()
{
    // Initial variables
    double rho_min = 0;
    double rho_max = 5.0;
    double epsilon = 1.0e-8;

    string outfilename = "Results_interacting_case.txt";
    ofile_interacting.open(outfilename);
    ofile_interacting << setiosflags(ios::showpoint | ios::uppercase);
    ofile_interacting << "Rho max is set to: " << rho_max << " and max error is set to: " <<
             epsilon << ". " << endl;
    ofile_interacting << " n:       Eigenvalues:                        " <<
             "Time (Armadillo):" << endl;

    // Loop over different n-values
    for(int i = 1; i < 2; i++)
    {
        int n = 2000 * i;
        mat A = zeros(n-1, n-1);
        //cout << "\nA(" << n-1 << "," << n-1 << "):\n" << A << endl;

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
        //cout << endl;
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
        //cout << "\nA(" << n-1 << "," << n-1 << "):\n" << A << endl;

        clock_t start_arma, finish_arma;
        start_arma = clock();
        cout << "start_arma: " << start_arma << " (in seconds: " <<
                (double)(start_arma/CLOCKS_PER_SEC) << ")" << endl;
        vec eigenval;
        mat eigenvec;

        mat B = A;
        eig_sym(eigenval, eigenvec, B);  // find eigenvalues/eigenvectors
        //cout << "\neigenval:\n" << eigenval << endl;
        eigenval = sort(eigenval);

        finish_arma = clock();
        cout << "finish_arma: " << finish_arma << " (in seconds: " <<
                (double)(finish_arma/CLOCKS_PER_SEC) << ")" << endl;
        double time_arma = (double)(finish_arma - start_arma) / CLOCKS_PER_SEC;
        cout << "time_arma: (s) " << time_arma << endl;

        ofile_interacting << setw(5) << n << setw(5);
        ofile_interacting << "{" << setprecision(8) << eigenval(0) << ", ";
        ofile_interacting << setprecision(8) << eigenval(1) << ", ";
        ofile_interacting << setprecision(8) << eigenval(2) << "}";
        ofile_interacting << setw(18) << setprecision(6) << time_arma << " s.";
        ofile_interacting << endl;
        cout << n << endl;
    }
    ofile_interacting.close();
}
