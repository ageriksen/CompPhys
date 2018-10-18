#include "earthsunsystem.h"

int main(int argc, char *argv[])
{
    mat r, v; 
    int t_end, N;
    double h, M;
    string method;
    t_end = 10; N = 3650; // timespan 10 years, 1 step per day
    h = (float)t_end/(float)N; //steplength
    M = 3e-6 // Msol, Earth mass in solar masses

    r = zeros(3, N);
    v = zeros(3, N);
    
    cout << "please choose a method: ";
    cin >> method;
    if( method == "VelocityVerlet" )
    {
        for(int i = 1; i < N; i++)
        {
            double ai = acceleration(r(:,i), M);
            r(i+1) = r(i) + h*v(i) + (h*h/2)*ai;
d:,ouble aii = acceleration(r(:,i+1));
            v(i+1) = v(i) = (h/2)*(aii + ai);
        }
    }
    else if( method == "ForwardEuler")
    { 
    }
    else
    {
        cout << "please ensure correct spelling." << endl;
        cout << "either 'VelocityVerlet' or 'ForwardEuler'" << endl;
        return 0;
    }
    return 0;
} //end of main

vec acceleration(vec r,double M)
{
    double R = norm(r);
    return 4*M_PI*r/(R*R*R);
}// end of acceleration
