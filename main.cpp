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

int main()
{
    jacobi_method();
    non_interacting();
    interacting();

    return 0;
}
