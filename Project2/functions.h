#ifndef functions_H
#define functions_H
#include <cmath>
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

void setup(int , double &, double & , double & , mat & , mat & );
void wrapper(double , mat & , mat & , int );
void Toeplitztridiag(mat &, int, double, double);
void jacobi_rotate( mat &, mat & , int & , int & , int  );
void offdiag( mat  , int & , int & , int );

void test_maxoffdiag();
void test_eigvalues();
void test_orthogonality();
#endif
