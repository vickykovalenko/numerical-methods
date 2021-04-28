#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

double it(double x)
{
    //return ((2*x-2)*x*x-36)/((3*x-4)*x-15);
    //return x - 0.5*(x*x*x - 2*x*x - 15*x + 36);
    //return 0.5*x*x*x-x*x-6.5*x+18;
    //return x + 0.38*(x*x*x + 3*x*x -1);
    //return x - 0.010596*(x*x*x+8*x*x-2);
    return x - 0.1052*(x*x*x+8*x*x-2);
}

double f(double x)
{
    //return ((x-2)*x-15)*x+36;
    return x*x*x+8*x*x-2;
}

int main()
{
    double x = 0, y = 1;
    while(fabs(x-y) > 0.0001)
    {
        cout<<x<<endl;
        y = x;
        x = it(x);
       
    }
    cout << "x = " << x << endl;
    cout << "f(x) = " << f(x) << endl;
}
