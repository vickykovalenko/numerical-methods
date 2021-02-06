//
//  main.cpp
//  сь
//
//  Created by Андрій on 2/5/21.
//  Copyright © 2021 Victoria. All rights reserved.
//

#include <iostream>
#include <cmath>
using namespace std;
double func(double x)
{
  return pow(x,5 )+ pow(x,2) + 3;
   
}

int main() {
    double a;
    double b;
    cin >> a;
    cin >> b;
    double mid;
    double eps = 0.0001;
    double x;
    while(b-a>eps)
    {
        if(func(a)*func(b)<0)
        {
            mid = (a+b)/2;
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
    x =( b+ a)/2;
    cout<<x<<endl;
    return 0;
}
