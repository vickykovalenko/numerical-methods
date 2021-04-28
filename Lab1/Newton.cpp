#include <stdio.h>
#include <math.h>
#include <iostream>

double fx(double x) { //return x*x-17;
    //return x*x-2;
    return sin(x);
}
double dfx(double x) { //return 2*x;
     //return 2*x;
    return cos(x);
}


double Newton(double x0)
{
    int count = 0;
    double x1  = x0 - fx(x0)/dfx(x0); //перше наближення
    
    while (fabs(x1-x0)>0.000001)
    {
        ++count;
        //std::cout<<x0<<std::endl;
        x0 = x1;
        x1 = x0 - fx(x0)/dfx(x0);
    }
    std::cout<<count<<std::endl;
    return x1;
}

int main () {
    double a = Newton(4);
    std::cout<<a<<std::endl;
    std::cout<<fx(a)<<std::endl;
    return 0;
}
    
    
    

