#include "earthsunsystem.h"

int main(int argc, char *argv[])
{
    vec r, v; 
    int t_end, N;
    double h, M;
    string method;
    t_end = 10; N = 3650; // timespan 10 years, 1 step per day
    h = (float)t_end/(float)N; //steplength

    r = zeros(
    
    cout << "please choose a method: ";
    cin >> method;
    if( method == "VelocityVerlet" )
    {
        for(int i = 1; i < N; i++)
        {
            double ai = Force(r(i))/M;
            r(i+1) = r(i) + h*v(i) + (h*h/2)*ai;
            double aii = Force(r(i+1));
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
