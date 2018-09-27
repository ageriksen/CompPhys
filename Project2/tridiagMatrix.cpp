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

void jacobi_rotate( mat & A, mat & R, int & k, int & l, int n )
{
    //declarations of sin and cos
    double s, c;
    // general jacobi_rotate algo
    if ( A(k,l) != 0.0 )
    {
        double t, tau; 
        tau = ( A(l,l) - A(k,k) )/( 2*A(k,l) );
        // avoiding loss of precision from difference between similar,
        // large numbers. 
        if ( tau >= 0 )
        {
            t = 1.0/( tau + sqrt( 1.0 + (tau*tau) ) );
        } 
        else
        {
            t = -1.0/( -tau + sqrt( 1.0 + (tau*tau) ) );
        }
        c = 1/sqrt( 1 + (t*t) );
        s = c*t;
    }
    else 
    {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = (c*c)*a_kk - 2.0*(c*s)*A(k,l) + (s*s)*a_ll; 
    A(l,l) = (s*s)*a_kk - 2.0*(c*s)*A(k,l) + (c*c)*a_ll; 
    // setting all nondiagonal elements manually
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    for ( int i = 0; i < n; i++)
    {
        if ( i != k && i != l)
        {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(i,l) = c*a_il - s*a_ik;
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);
        }
        // adjusting the eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);

        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
} //end of jacobi_rotate

// off diagonal function using armadillo to find largest 
// offdiagonal element
void offdiag( mat & A, int & p, int & q, int n )
{
    double max; 
    for ( int i = 0; i < n; i++)
    {
        for ( int j = i+1; j < n; j++)
        {
            double aij = fabs( A(i,j) );
            if ( aij > max )
            {
                max = aij; p = i; q = j;
            }
        }
    }
}
