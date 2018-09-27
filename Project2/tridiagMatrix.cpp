#include <cmath>
#include <armadillo>

using namespace arma;
using namespace std;

//function declarations

int main()
{
    //initializing variables
    int n = 3; // dim of matrix. chosen to ensure validity before 
               // generalizing to larger matrices
    double h = (double)(1.0/n); 
    double d = 2.0/(h*h); 
    double a = -1.0/(h*h);
    double eps_avg;
   
    //creating and filling matrices
    mat A = zeros(n, n); //what will be our tridiagonal Toeplitz matrix
    vec lmbd = zeros(n); //our analytical eigenvalues
    vec eps = zeros(n); //relative error between numerical and 
                        // analytical eigenvalues. 
    for( int i = 1; i < (n-1); i++)
    {
        // filling Toeplitz matrix
        A(i,i) = d; A(i, i+1) = a; A(i, i-1) = a;
        // picking analytical eigenvalues
        lmbd(i) = d + 2*a*cos(i*M_PI/(n+1));
    }

    // supplying the ends with appropriate appropriate  elements
    A(0,0) = d; A(0,1) = a; A(n-1,n-1) = d; A(n-1,n-2) = a;

}

void jacobi_method(double **A, double **R, int n)
{
    //set up eig vector matrix
    for ( int i = 0; i<n; i++){
        for ( int j = 0; i<n; j++){
            if ( i == j){
                R[i][j] = 1.0; 
            } else {
                R[i][j] = 0.0;
            }
        }
    }

    int k, l;
    double epsilon = 1.0e-8;
    double max_number_iterations = (double) n * (double) n * (double) n;
    int iteration = 0;
    double max_offdiag = maxoffdiag(A, &k, &l, n);
    
    while (fabs(max_offdiag) > epsilon 
            && (double) iteration < max_numeber_iterations){
}

// finding the maximum offdiagonal matrix element. optimizable?
double maxoffdiag( double **A, int *k, int *l, int n)
{
    double max = 0.0;
    
    for ( int i = 0; i < n; i++){
        for ( int j = i+1; j < n; i++){
           if ( fabs(A[i][j]) > max){
              max = fabs(A[i][j]);
              *l = i;
              *k = k;
           }
        }
    }
    return max;
}

}
