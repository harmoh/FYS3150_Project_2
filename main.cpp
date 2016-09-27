#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <string>
#include <time.h>
#include "jacobi.h"
#include "interacting_case.h"

using namespace std;

int main()
{
    //jacobi_method();

    clock_t start = clock();
    interacting();
    clock_t end = clock();
    double total_time = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "Total time (s): " << total_time << endl;

    return 0;
}
