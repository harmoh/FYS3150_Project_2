#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <string>
#include <time.h>
#include "jacobi.h"
#include "non_interacting_case.h"
#include "interacting_case.h"

using namespace std;
using namespace arma;

int main()
{
    //jacobi_method();

    // Unit test
    // Create and initialze a nxn test matrix
    int n = 6;
    int k = 0;
    int l = 0;
    mat A = zeros(n-1, n-1);
    A(n-2, n-2) = 10;
    cout << A << endl;
    for(int i = 0; i < n-2; i++)
    {
        A(i,i) = 2 * i + 10;
        A(i,i+1) = i * (1 - i) + 3;
        A(i+1,i) = i * i + 5;
    }
    double max = max_offdiagonal(n-1, A, &k, &l);

    cout << A << endl << "Max: " << max << endl;

    //non_interacting();
    //interacting();

    return 0;
}
