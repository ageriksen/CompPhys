#include "functions.h"

int main(int argc, char *argv[])
{
    if ( argc < 2 )
    {
        cout    << "please supply a dimension for matrix size and integration step length."
                << "as well as the name of file written, as 'n'DimEigval.txt"
                << endl; 
        return 1;
    }

    // initializing variables
    int maxdistance, iterations, dim; 
    double diagonal, semidiagonal, tolerance, step, time, sum_off, eps_rel1, eps_rel2, eps_rel3;
    mat eigvalmatrix, eigvecmatrix;
    vec outvalues, eigenvalues, eigenstate1, eigenstate2, eigenstate3, omega; 
    uvec eigval_index; 
    
    // establishing necessary elements
    maxdistance = 10; // maxdistance rho = r/alpha
    diagonal = 2.0; semidiagonal = -1.0;
    tolerance = 1.0e-10; dim = atoi(argv[1]);
    iterations = 0; omega = zeros(4);
    omega(0) = 0.01; omega(1) = 0.5; omega(2) = 1.0; omega(3) = 5;

    for ( unsigned i = 0; i < omega.size(); i++)
    {
        // setup and execution of jacobi method
        setup(dim, maxdistance,  step, diagonal, semidiagonal, eigvalmatrix, eigvecmatrix, omega(i));
        wrapper( tolerance, iterations, time, eigvalmatrix, eigvecmatrix, dim);
        sum_offdiag(eigvalmatrix, dim, sum_off);
        
        // finding the lowest eigenstates of the HO system.
        eigenvalues = zeros(dim);
        for ( unsigned i = 0; i < eigenvalues.size(); i++)
        {
            eigenvalues(i) = eigvalmatrix(i, i);
        }
        eigval_index = sort_index(eigenvalues);
        eigenstate1 = eigvecmatrix.col(eigval_index(0));
        eigenstate2 = eigvecmatrix.col(eigval_index(1));
        eigenstate3 = eigvecmatrix.col(eigval_index(2));

        // investigating whether or not the lowest 3 eigenvalues are correct to within
        // 4 leading digits. 
        eps_rel1 = ( eigenvalues(eigval_index(0)) - 3.0 )/3.0;
        eps_rel2 = ( eigenvalues(eigval_index(1)) - 7.0 )/7.0;
        eps_rel3 = ( eigenvalues(eigval_index(2)) - 11.0 )/11.0;
        // reporting error
        cout << "relative error for orbital states 1 2 and 3:" << endl;
        cout << eps_rel1 << setw(18) << eps_rel2 << setw(18) << eps_rel3 << endl;

        //print results
        cout    << "run of HO potential omega"<<i+1<<" = "<< omega(i) << "for \n"
                << "A rotation of a matrix with dimensions n: " << dim
                << " diagonalizes in " << iterations << "rotations.\n"
                << " These rotations are completed over a period of " << time << " seconds"
                << endl;

        // write result to file:
        string outfile = "omega"+to_string(i)+"dim"+(string)argv[1]+"Eigval.txt";
        ofstream myfile;
        myfile.open(outfile, ios::out);
        myfile << "eigenstate1" << setw(20)
               << "eigenstate2" << setw(20)
               << "eigenstate3" << setw(20)
               << "rotations" << setw(20) 
               << "time" << setw(20) 
               << "ofdiagonal sum" << setw(20)
               << "number of steps" << endl;
        myfile << eigenstate1(0) << setw(20)
               << eigenstate2(0) << setw(20)
               << eigenstate3(0) << setw(20)
               << iterations << setw(20) 
               << time << setw(20) 
               << sum_off << setw(20)
               << dim << endl;
        for ( unsigned i = 1; i < eigenstate1.size(); i++)
        {
            myfile << eigenstate1(i) << setw(20)
                   << eigenstate2(i) << setw(20)
                   << eigenstate3(i) << endl;
        }
        myfile.close();
    }
    // if test to check if tests should be run.
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
        else if ( argument == 3 )
        {
            cout << "3 is a test to confirm 1 rotation for a 2 dim matrix" << endl;
            test_2dimrotation();
        }
        else if ( argument == 4 )
        {
            cout << "4 to test rotations to matrix dimensions" << endl;
            test_jacobiiterations();
        }
     }
     return 0;
} // end of main function

//############################################################
//setup function to set up necessary variables.
void setup(int dim, int maxdistance, double & step, double & diagonal, double & semidiagonal, mat & eigvalmatrix, mat & eigvecmatrix, double omega)
{
    //setting diagonal and semidiagonal element values
    double diag, semidiag; 
    diag = diagonal; semidiag = semidiagonal;
    step = ((double)maxdistance/dim); // step size
    diagonal = diag/(step*step);  // diagonal constant from taylor expand of diff-eq.
    semidiagonal = semidiag/(step*step); //  constant for future and past step wit taylor expansion
    //filling matrices
    eigvalmatrix = zeros<mat>(dim, dim); //what will be our tridiagonal Toeplitz matrix
    eigvecmatrix = eye<mat>(dim,dim); // identity matrix to contain our eigenvectors
    Toeplitztridiag( eigvalmatrix, dim, step, diagonal, semidiagonal, omega); // tridiagonalizing eigenvalues
} // end of setup

//#############################################################
//Wrapper function to loop over iterations of rotations
void wrapper(double tolerance,int & iterations, double & time,  mat & eigvalmatrix, mat & eigvecmatrix, int dim )
{
    double maxnondiagonal = 1; 
    auto start = chrono::system_clock::now();
    int p, q; 
    while ( maxnondiagonal > tolerance )
    {
        offdiag(eigvalmatrix, p, q, dim);
        maxnondiagonal = fabs(eigvalmatrix(p,q));
        jacobi_rotate( eigvalmatrix , eigvecmatrix , p, q, dim);
        iterations ++;
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> runtime = end - start;
    time = runtime.count();
}

//#############################################################
//function to set a tridiagonal Toeplitx matrix given diagonal 
//and "neighbouring" elements values d and a
void Toeplitztridiag(mat & Matrix, int dim, double step, double diagonal, double semidiagonal, double omega)
{
    double rho;
    for( int i = 1; i < (dim-1); i++)
    {
        // filling Toeplitz matrix((i+1)*step)*((i+1)*step); 
        rho = (double)i*step;
        Matrix(i,i) = diagonal + omega*(rho*rho) + 1.0/rho;
        Matrix(i, i+1) = semidiagonal; 
        Matrix(i, i-1) = semidiagonal;
    }

    // supplying the ends with appropriate appropriate  elements
    rho = step;
    Matrix(0,0) = diagonal + omega*(rho*rho) + 1.0/rho;  
    Matrix(0,1) = semidiagonal; 
    rho = (double)(dim-1)*step;
    Matrix(dim-1,dim-1) = diagonal + omega*(rho*rho) + 1.0/rho; 
    Matrix(dim-1,dim-2) = semidiagonal;
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
    double aij; 
    for ( int i = 0; i < n; i++)
    {
        for ( int j = i+1; j < n; j++)
        {
            aij = fabs( A(i,j) );
            if ( aij > max )
            {
                max = aij; p = i; q = j;
            }
        }
    }
} // end of offdiag

void sum_offdiag( mat eigenvalues, int dim, double & sum_off)
{
    for ( int i = 0; i < dim; i++)
    {
        for ( int j = i+1; j < dim; j++)
        {
            sum_off += eigenvalues(i, j);
        }
    }

} // end of sum_offdiag(...)

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
    int dim, maxdistance, iterations;
    double diagonal, semidiagonal, tolerance, epsilon, step, time;
    mat testeigvalmatrix, testeigvecmatrix;
    vec lambda, testeigvalvector;
    dim = 5; maxdistance = 1;  diagonal = 2.0; semidiagonal = -1.0; 
    tolerance = 1.0e-15; iterations = 0;
    
    setup(dim, maxdistance, step, diagonal, semidiagonal, testeigvalmatrix, testeigvecmatrix, 1);

    wrapper( tolerance, iterations, time,  testeigvalmatrix, testeigvecmatrix, dim);
    
    lambda = zeros<vec>(dim);
    testeigvalvector = zeros<vec>(dim);
    for ( int i = 0; i < dim; i++)
    {
        lambda(i) = diagonal + 2*semidiagonal*cos( (double(i)+1)*M_PI /( double(dim)+1 ));
        testeigvalvector(i) = testeigvalmatrix(i,i);
    }
    cout << lambda << endl;
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
    int dim, maxdistance, iterations;
    double diagonal, semidiagonal, tolerance, testinnerproduct, step, time;
    mat testeigvalmatrix, testeigvecmatrix;
    vec testeigveccollumn_i, testeigveccollumn_j;
    dim = 5; maxdistance = 1; diagonal = 2.0; semidiagonal = -1.0; 
    tolerance = 1.0e-15; iterations = 0;

    setup( dim, maxdistance, step, diagonal, semidiagonal, testeigvalmatrix, testeigvecmatrix, 1);

    testeigveccollumn_i = testeigvecmatrix.col(1);
    testeigveccollumn_j = testeigvecmatrix.col(3);
    testinnerproduct = dot(testeigveccollumn_i, testeigveccollumn_j);

    cout << "our initial Toeplitz matrix has an innerproduct between two collumns of" << endl;
    cout << testinnerproduct << endl;

    wrapper( tolerance, iterations, time,  testeigvalmatrix, testeigvecmatrix, dim );

    testeigveccollumn_i = testeigvecmatrix.col(1);
    testeigveccollumn_j = testeigvecmatrix.col(3);
    testinnerproduct = dot(testeigveccollumn_i, testeigveccollumn_j);

    cout << "our rotated matrix has an innerproduct between two collumns of" << endl;
    cout << testinnerproduct << endl;
} // end of test_orthogonality

// ############ TEST IDEA: CHECK THAT 2X2 MATRIX NEEDS ONLY 1 ROTATION TO DIAGONALIZE
// void test_2dimrotation()
void test_2dimrotation()
{
    int iterations, dim;
    double diagonal, semidiagonal, tolerance, step, time;
    mat eigvalmatrix, eigvecmatrix; 
    diagonal = 2.0; semidiagonal = -1.0;
    tolerance = 1.0e-16; dim = 2; iterations = 0;
    setup(dim, 1, step, diagonal, semidiagonal, eigvalmatrix, eigvecmatrix, 1);
    wrapper( tolerance, iterations, time, eigvalmatrix, eigvecmatrix, dim);
    cout << "test of 2 dim rotation reveals " << iterations 
         << " rotations before we have the eigenvalues" << endl;
} // end of test_2dimrotation()    

void test_jacobiiterations()
{
    int maxdistance, iterations, dim;
    double diagonal, semidiagonal, tolerance, step, time;
    mat eigvalmatrix, eigvecmatrix;
    vec outvalues;

    maxdistance = 0; // maxdistance rho = r/alpha
    diagonal = 2.0; semidiagonal = -1.0;
    tolerance = 1.0e-5; 

    vec repeatarray = linspace(2, 100, 99);
    vec rotations = zeros(repeatarray.size());
    vec timed = zeros(repeatarray.size());
    for (unsigned i = 0; i < (repeatarray.size() - 1); i++ )
    {
        dim = repeatarray(i);
        iterations = 0; 
        setup(dim, maxdistance,  step, diagonal, semidiagonal, eigvalmatrix, eigvecmatrix, 1);
        wrapper( tolerance, iterations, time, eigvalmatrix, eigvecmatrix, dim);
        rotations(i) = iterations;
        timed(i) = time;
    }
    ofstream myfile;
    myfile.open("DimRotationTime.txt");
    myfile << "dimensions" << setw(20) 
           << "rotations" << setw(20)
           << "time spent" << endl;
    for ( unsigned i = 0; i < rotations.size(); i++ )
    {
        myfile << repeatarray(i) << setw(20) << rotations(i) << setw(20) << timed(i) << endl;
    }
    myfile.close();
} // end of test_jacobiiterations
