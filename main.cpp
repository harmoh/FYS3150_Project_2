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

    for(int i = 1; i < 8; i++)
    {
        int n = 50 * i;
        jacobi_method(n);
    }

    return 0;
}
