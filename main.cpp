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
    //jacobi_method();

    test_max_off();
    test_orthogonal();

    non_interacting();
    interacting();

    return 0;
}
