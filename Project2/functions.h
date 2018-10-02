#ifndef functions_H
#define functions_H
#include <cmath>
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace arma;

void setup(int , int , double &, double & , double & , mat & , mat & );
void wrapper(double , int &, double & , mat & , mat & , int);
void Toeplitztridiag(mat &, int, double,  double, double);
void jacobi_rotate( mat &, mat & , int & , int & , int  );
void offdiag( mat  , int & , int & , int );
void sum_offdiag( mat , int , double &);

void test_maxoffdiag();
void test_eigvalues();
void test_orthogonality();
void test_2dimrotation();
void test_jacobiiterations();
#endif
