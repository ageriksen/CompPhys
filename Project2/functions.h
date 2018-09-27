#ifndef functions_H
#define functions_H
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

void Toeplitztridiag(mat &, int, double, double);
void jacobi_rotate( mat &, mat & , int & , int & , int  );
void offdiag( mat  , int & , int & , int );
#endif
