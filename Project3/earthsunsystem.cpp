#include "earthsunsystem.h"

int main(int argc, char *argv[])
{
    mat r, v; 
    int t_end, N, method;
    double h, M;
    t_end = 10; N = 3650; // timespan 10 years, 1 step per day
    h = (float)t_end/(float)N; //steplength
    M = 3e-6; // Msol, Earth mass in solar masses

    r = zeros(3, N);
    v = zeros(3, N);
    r(0,0) = 1;
    v(1,0) = 2*M_PI; // AU/yr velocity of earth. 1 period = 1 year, r = 1AU. 
    
    
    cout << "please choose a method, either 'VelocityVerlet'[0] or 'ForwardEuler'[1]: ";
    cin >> method;
    if( method == 0)
    {
        for(int i = 1; i < N-1; i++)
        {
            vec ai = acceleration(r.col(i), M);
            r.col(i+1) = r(i) + h*v.col(i) + (h*h/2)*ai;
            vec aii = acceleration(r.col(i+1), M);
            v.col(i+1) = v.col(i) + (h/2)*(aii + ai);
        }
    }
    else if( method == 1 )  
    { 
        for(int i = 1; i < N-1; i++)
        {
            v.col(i+1) = v.col(i) + acceleration(r.col(i), M)*h;
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
    cout << "do you wish to write results(y, n)? ";
    cin >> reply;
    if( reply == "y")
    {
        string filename;
        cout << "please provide a filname: ";
        cin >> filename;
        ofstream myfile;
        myfile.open(outfile, ios::out);
        myfile << "x" << setw(20) << "y" << setw(20) << "z" << endl;
        for ( unsigned i = 1; i < )
    }


    return 0;
} //end of main

vec acceleration(vec r,double M)
{
    double R = norm(r);
    return 4*M_PI*r/(R*R*R);
}// end of acceleration
