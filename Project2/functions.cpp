#include "functions.h"

int main(int argc, char *argv[])
{
    if ( argc < 2 )
    {
        cout    << "please supply a filename for storage "
                << "and optionally numbers for corresponding unit tests"
                << endl; 
        return 1;
    }
    int dim = 5;
    double diagonal, semidiagonal, tolerance, step;
    mat eigvalmatrix, eigvecmatrix;
    diagonal = 2.0; semidiagonal = -1.0;
    tolerance = 1.0e-10;
 
    setup(dim, step, diagonal, semidiagonal, eigvalmatrix, eigvecmatrix);
    wrapper( tolerance, eigvalmatrix, eigvecmatrix, dim);
    
    if (argc > 2)
    {
        int argument = atoi(argv[2]);
        cout << "additional raguments mean initiating tests. " << endl;
        if (argument == 0)
        {
            cout << "0 tests that max values are properly picked." << endl;
            test_maxoffdiag();
        }
        else if (argument == 1)
        {
            cout << "1 tests conservation of eigenvalues under jacobi transformations" <<
                "of Toeplitz matrix." << endl;
            test_eigvalues();
        }
        else if (argument == 2)
        {
            cout << "2 tests orthogonality of eigenvectors under rotation." << endl;
            test_orthogonality();
        }
     }
     return 0;
} // end of main function

//############################################################
//setup function to set up necessary variables.
void setup(int dim, double & step, double & diagonal, double & semidiagonal, mat & eigvalmatrix, mat & eigvecmatrix )
{
    //setting diagonal and semidiagonal element values
    double diag, semidiag; 
    diag = diagonal; semidiag = semidiagonal;
    step = (double)(1.0/dim); // step size
    diagonal = diag/(step*step);  // diagonal constant from taylor expand of diff-eq.
    semidiagonal = semidiag/(step*step); //  constant for future and past step wit taylor expansion
    //filling matrices
    eigvalmatrix = zeros<mat>(dim, dim); //what will be our tridiagonal Toeplitz matrix
    eigvecmatrix = eye<mat>(dim,dim); // identity matrix to contain our eigenvectors
    Toeplitztridiag( eigvalmatrix, dim, step, diagonal, semidiagonal ); // tridiagonalizing eigenvalues
} // end of setup

//#############################################################
//Wrapper function to loop over iterations of rotations
void wrapper(double tolerance, mat & eigvalmatrix, mat & eigvecmatrix, int dim )
{
    double maxnondiagonal = 1; 
    while ( maxnondiagonal > tolerance )
    {
        int p, q; 
        offdiag(eigvalmatrix, p, q, dim);
        maxnondiagonal = fabs(eigvalmatrix(p,q));
        jacobi_rotate( eigvalmatrix , eigvecmatrix , p, q, dim);
    }
}

//#############################################################
//function to set a tridiagonal Toeplitx matrix given diagonal 
//and "neighbouring" elements values d and a
void Toeplitztridiag(mat & Matrix, int dim, double step, double diagonal, double semidiagonal)
{
    for( int i = 1; i < (dim-1); i++)
    {
        // filling Toeplitz matrix
        Matrix(i,i) = diagonal + (i*step)*(i*step); 
        Matrix(i, i+1) = semidiagonal; 
        Matrix(i, i-1) = semidiagonal;
    }

    // supplying the ends with appropriate appropriate  elements
    Matrix(0,0) = diagonal; Matrix(0,1) = semidiagonal; 
    Matrix(dim-1,dim-1) = diagonal + ((dim-1)*step)*((dim-1)*step); Matrix(dim-1,dim-2) = semidiagonal;
} // end of Toeplitztridiag

//#############################################################
//function using orthogonal rotations of collumns to produce a
//tridiagonal matrix within a certain tolerance.
void jacobi_rotate( mat & A, mat & R, int & k, int & l, int n )
{
    //declarations of sin and cos
    double s, c;
    
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
    A(l,l) = (s*s)*a_kk + 2.0*(c*s)*A(k,l) + (c*c)*a_ll; 
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
            A(i,l) = c*a_il + s*a_ik;
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

//####################################################
// function using armadillo to find largest 
// offdiagonal element
void offdiag( mat A, int & p, int & q, int n )
{
    double max = 0.0; 
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
} // end of offdiag

// ########################################################
// ########################################################
// TESTS

// ################################
// test of offdiagonal max value function
void test_maxoffdiag()
{
   /*
    *Function to test that the "offdiag" function actually picks 
    out the maximum value of the matrix sent in. To do this, I 
    send in a known 5x5 matrix and ensure that the maximum value
    picked is the correct one. 
    */
    int test_p, test_q, dim, maxvalue;
    dim = 5;
    maxvalue = 10;
    mat testmatrix;
    testmatrix = eye<mat>(dim,dim);
    testmatrix(dim - 3, dim -2) = maxvalue;
    offdiag( testmatrix, test_p, test_q, dim);
    cout << "test of 'offdiag' function using "<<
       "simple matrix of dim " << dim << endl;
    cout << "I chose a maxvalue of "<< maxvalue << 
        " and assigned it to testmatrix("<< dim -3 << ","<< dim -2 << ")" << endl;
    cout << "offdiag returns values of p, q as " << test_p << ", " << test_q << endl;
    cout << "which corresponds to a maxvalue of " << testmatrix(test_p, test_q) << endl;
} // end of test_maxoffdiag

// ##################################
// test that correct eigenvalues extracted. 
void test_eigvalues()
{
    /*
     * This test takes a jacobi rotated matrix and compares the 
     * eigenvalues along the diagonal and compares it with the 
     * analytical eigenvalues from the fomula bellow eq (2) in 
     * the project 2 pdf. 
     */
    int dim;
    double diagonal, semidiagonal, tolerance, epsilon, step;
    mat testeigvalmatrix, testeigvecmatrix;
    vec lambda, testeigvalvector;
    dim = 5; diagonal = 2.0; semidiagonal = -1.0; 
    tolerance = 1.0e-15; 
    
    setup(dim, step, diagonal, semidiagonal, testeigvalmatrix, testeigvecmatrix);

    wrapper( tolerance, testeigvalmatrix, testeigvecmatrix, dim);
    
    lambda = zeros<vec>(dim);
    testeigvalvector = zeros<vec>(dim);
    for ( int i = 0; i < dim; i++)
    {
        lambda(i) = diagonal + 2*semidiagonal*cos( (double(i)+1)*M_PI /( double(dim)+1 ));
        testeigvalvector(i) = testeigvalmatrix(i,i);
    }
    epsilon = (mean(testeigvalvector)-mean(lambda))/mean(lambda); // relative difference between means
                                                                  // of analyrical and numerical eigvals.
    
    cout << "test eigenvalue matrix" << endl;
    cout << testeigvalmatrix << endl;
    cout << "analytical eigenvalue vector lambda: " << endl;
    cout << lambda << endl;
    cout << "the average error between analytical and numerical eigenvalues are:" <<
        epsilon << endl;
} // end of test_eigvalues

// ###########################################
// test that jacobi rotations maintain orthogonality.
void test_orthogonality()
{
    /*
     * The test aims to ensure that each column in the transformed matrix maintains
     * orthogonality under rotations. The test takes a small matrix and ensures that 
     * columns are orthogonal before and after rotation.
     */
    int dim;
    double diagonal, semidiagonal, tolerance, testinnerproduct, step;
    mat testeigvalmatrix, testeigvecmatrix;
    vec testeigveccollumn_i, testeigveccollumn_j;
    dim = 5; diagonal = 2.0; semidiagonal = -1.0; 
    tolerance = 1.0e-15; 

    setup( dim, step, diagonal, semidiagonal, testeigvalmatrix, testeigvecmatrix);

    testeigveccollumn_i = testeigvecmatrix.col(1);
    testeigveccollumn_j = testeigvecmatrix.col(3);
    testinnerproduct = dot(testeigveccollumn_i, testeigveccollumn_j);

    cout << "our initial Toeplitz matrix has an innerproduct between two collumns of" << endl;
    cout << testinnerproduct << endl;

    wrapper( tolerance, testeigvalmatrix, testeigvecmatrix, dim );

    testeigveccollumn_i = testeigvecmatrix.col(1);
    testeigveccollumn_j = testeigvecmatrix.col(3);
    testinnerproduct = dot(testeigveccollumn_i, testeigveccollumn_j);

    cout << "our rotated matrix has an innerproduct between two collumns of" << endl;
    cout << testinnerproduct << endl;
} // end of test_orthogonality

// ############ TEST IDEA: CHECK THAT 2X2 MATRIX NEEDS ONLY 1 ROTATION TO DIAGONALIZE
// void test_2dimrotation()
// {
//      int dim = 2;
//      double ...;
//      mat ...; 
//      iterations = 0;
//      setup(...); 
//      wrapper(...); 
//      if iteration > 1:
//          break, print( "too many rotations to find eigenvals for 2x2 matrix)
// } // end of test_2dimrotation
