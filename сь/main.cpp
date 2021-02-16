#include <iostream>
#include <cmath>
using namespace std;
double func(double x)
{
  //return pow(x,5)+ pow(x,2) + 3;
    return x*x*x - x*x + 2;

   
}

int main() {
    double a;
    double b;
    cin >> a;
    cin >> b;
    double mid = 0;
    double eps = 0.00001;
    //double x;
    if(func(a)==0)
    {
        cout<<"root:"<<a<<endl;
        return 0;
    }
    if(func(b)==0)
    {
        cout<<"root:"<<b<<endl;
        return 0;
    }
    while((b-a>eps)||func(mid)==0)
    {
        mid = (a+b)/2;

        if(func(a)*func(b)<0)
        {
            if(func(a)*func(mid)<0)
            {
                b = mid;
            }
            else
            {
                a = mid;
            }
        }
        else {
            cout<<"no roots here"<<endl;
            return 0;
        }
    }
    //x =( b+ a)/2;
    cout<<"root "<<mid<<endl;
    return 0;
}
