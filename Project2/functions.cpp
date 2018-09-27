#include "functions.h"


//function declarations

int main()
{
    //initializing variables
    int n = 3; // dim of matrix. chosen to ensure validity before 
               // generalizing to larger matrices
    double h, d, a; 
    h = (double)(1.0/n); 
    d = 2.0/(h*h); 
    a = -1.0/(h*h);
    //creating and filling matrices
    mat A, R;
    A = zeros(n, n); //what will be our tridiagonal Toeplitz matrix
    R = zeros(n, n); //our eigval matrix
    Toeplitztridiag( A, n, d, a );
    Toeplitztridiag( R, n, d, a );

    double tolerance = 1.0e-10; 
    double maxnondiagonal = tolerance*10;
    int iterations = 0;
    while ( maxnondiagonal > tolerance && iterations < n ) 
    {
        int p, q; 
        offdiag(A, p, q, n);
        jacobi_rotate( A , R , p, q, n);
        maxnondiagonal = A(p,q);
        iterations++;
    }

    return 0;
}

void Toeplitztridiag(mat & Matrix, int n, double d, double a)
{
    for( int i = 1; i < (n-1); i++)
    {
        // filling Toeplitz matrix
        Matrix(i,i) = d; Matrix(i, i+1) = a; Matrix(i, i-1) = a;
    }

    // supplying the ends with appropriate appropriate  elements
    Matrix(0,0) = d; Matrix(0,1) = a; 
    Matrix(n-1,n-1) = d; Matrix(n-1,n-2) = a;
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
