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
#include "unit_test.h"

using namespace std;
using namespace arma;

int main()
{
    cout << "Performing Jacobi rotation..." << endl;
    jacobi_method();

    cout << "Unit testing..." << endl;
    test();

    cout << "Non-interacting and interacting...";
    non_interacting();
    interacting();
    cout << "Finished!" << endl;

    return 0;
}
