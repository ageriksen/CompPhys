#include "JacobiMethodEigenvalues.h"

int main()
{
    cout << "hello world!" << endl;
} // end of main(...)

// #### jacobi rotation algorithm ####
void Rotate( mat & eigenvalues, mat & eigenvectors, int & k, int & l, int dim)
{
    //setup of sin and cosine for orthogonal rotations
    double s, c;
    if ( eigenvalues(k, l) != 0.0)
    {        
        // initiate tangent, t, and Tau
        double t, tau;
        tau = ( eigenvalues(l, l) - eigenvalues(k, k) )/( 2*eigenvalues(k, l) );
        // to avoid loss of precision at difference of large nr, 
        // t = -tau +/- sqrt( 1 + tau^2 ) *(t_(c.c)/t(c.c))
        if ( tau >= 0)
        {
            t = 1.0/( tau + sqrt( 1.0 + (tau*tau) ) );
        }
        else
        {
            t = - 1.0/( -tau + sqrt( 1.0 + (tau*tau) ) );
        }
        // picking c and s from relation to tangent t
        c = 1.0/sqrt( 1.0 + (t*t) );
        s = c*t;
    }
    else
    {
        c = 1.0;
        s = 0.0;
    }

    // rotating
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = eigenvalues(k,k); a_ll = eigenvalues(l,l);
    eigenvalues(k,k) = (c*c)*a_kk - 2.0*(c*s)*eigenvalues(k,l) + (s*s)*a_ll;
    eigenvalues(l,l) = (s*s)*a_kk + 2.0*(c*s)*eigenvalues(k,l) + (c*c)*a_ll;
    // forcing nondiagonals 0
    eigenvalues(k,l) = 0.0; eigenvalues(l,k) = 0.0;
    for ( int i = 0; i < dim; i++ )
    {
        if ( i != k && i != l )
        {
            a_ik = eigenvalues(i,k); a_il = eigenvalues(i,l);
            eigenvalues(i,k) = c*a_ik - s*a_il;
            eigenvalues(i,l) = c*a_il + s*a_ik;
            eigenvalues(k,i) = eigenvalues(i,k); 
            eigenvalues(l,i) = eigenvalues(i,l);
        }
        //adjusting eigevectors
        r_ik = eigenvectors(i,k); r_il = eigenvectors(i,l);
        eigenvectors(i,k) = c*r_ik - s*r_il;
        eigenvectors(i,l) = c*r_il + s*r_ik;
    }
} // end of Rotate(...)
