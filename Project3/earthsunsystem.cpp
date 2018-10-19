#include "earthsunsystem.h"

int main(int argc, char *argv[])
{
    mat r, v; 
    int t_end, N, method;
    double h;
    t_end = 10; N = 1000; // timespan in years and number of timesteps
    h = (float)t_end/(float)N; //steplength

    r = zeros(2, N);
    v = zeros(2, N);
    r(0,0) = 1;
    v(1,0) = 2*M_PI; // AU/yr velocity of earth. 1 period = 1 year, r = 1AU. 
    
    
    cout << "please choose a method, either 'VelocityVerlet'[0] or 'ForwardEuler'[1]: ";
    cin >> method;
    if( method == 0) // Velocity Verlet
    {
        for(int i = 0; i < N-1; i++)
        {
            vec ai = acceleration(r.col(i));
            r.col(i+1) = r.col(i) + h*v.col(i) + h*h*ai/2;
            vec aii = acceleration(r.col(i+1));
            v.col(i+1) = v.col(i) + (h/2)*(aii + ai);
        }
    }
    else if( method == 1 )   //Forward Euler
    { 
        for(int i = 0; i < N-1; i++)
        {
            v.col(i+1) = v.col(i) + acceleration(r.col(i))*h;
            r.col(i+1) = r.col(i) + v.col(i+1)*h;
        }
    }
    else
    {
        cout << "please ensure correct spelling." << endl;
        cout << "either 'VelocityVerlet' or 'ForwardEuler'" << endl;
        return 0;
    }

    //writing
    string reply;
    int skip;
    cout << "do you wish to write results(y, n)? ";
    cin >> reply;
    if( reply == "y")
    {
        string filename;
        cout << "please provide a filname with extension: ";
        cin >> filename;
        cout << "what increment do you wish to write? ";
        cin >> skip;
        ofstream myfile;
        myfile.open(filename, ios::out);
        myfile << "#x" << setw(20) << "y \n";
        for ( unsigned i = 0; i < r.n_cols; i = i + skip)
        {
            myfile <<
                r(0,i) << setw(20) << r(1,i) << setw(20) << "\n";
        }
        myfile.close();
    }


    return 0;
} //end of main

vec acceleration(vec r)
{
    double R = norm(r);
    return -4*M_PI*M_PI*r/(R*R*R);
}// end of acceleration

void verlet

void euler 
