#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <string>
#include "jacobi.h"

using namespace std;
using namespace arma;

int main()
{
    int n = 10;

    jacobi_method(n);

    return 0;
}
