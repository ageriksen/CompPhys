#include <cmath>
#include <armadillo>

using namespace arma;
//using namespace std;

int main()
{
    //initializing variables
    int n = 10;
    double h = (double)(1.0/n);
    double d = 2.0/h;
    double a = -1.0/h;
   
    //creating and filling matrices
    mat A = zeros(n, n);
    vec lmbd = zeros(n);
    vec eps = zeros(n);
    for( int i = 1; i < (n-1); i++){
        A(i,i) = d;
        A(i, i+1) = a;
        A(i, i-1) = a;
        // picking analytical eigenvalues
        lmbd(i) = d + 2*a*cos(i*M_PI/(n+1));
    }

    // supplying remaining tridiagonal elements
    A(0,0) = d;
    A(0,1) = a;
    A(n-1,n-1) = d;
    A(n-1,n-2) = a;
    
    // extracting numerical eigenvalues with armadillo
    vec eigval = eig_sym(A);
    eigval.print();

    for( int i = 1; i < n; i++){
        eps(i) = fabs(eigval(i) - lmbd(i)) / lmbd(i);
    }
}
