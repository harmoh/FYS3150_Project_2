#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <string>
#include "jacobi.h"

using namespace std;

int main()
{
    int n = 50;

    jacobi_method(n);

    return 0;
}
