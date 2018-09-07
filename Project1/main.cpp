#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    mat A = randu<mat>(4,5);
    cout << A << endl;
    return 0;
}
