#ifndef JACOBI_METHOD_H
#define JACOBI_METHOD_H

void jacobi_method();
double max_offdiagonal(int n, arma::mat A, int *k, int *l);
void rotate(int n, arma::mat &A, arma::mat &R, int k, int l);
void eigenvalues_arma(arma::mat &A);

#endif // JACOBI_METHOD_H
