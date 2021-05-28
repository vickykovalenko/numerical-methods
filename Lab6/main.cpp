#include <iostream>
#include <math.h>
using namespace std;

#define eps 0.00001

double function1(double x, double y)
{
    return pow(x, 2) / pow(y, 2) - cos(y) - 2;
}

double function2(double x, double y)
{
    return pow(x, 2) + pow(y, 2) - 6;
}

double func1dx(double x, double y)
{
    return 2 * x / pow(y, 2);
}

double func1dy(double x, double y)
{
    return -2 * pow(x, 2) / pow(y, 3) + sin(y);
}

double func2dx(double x, double y)
{
    return 2 * x;
}

double func2dy(double x, double y)
{
    return 2 * y;
}

double function3(double x, double y)
{
    return x * x * x + y * y - 1 - 4;
}
double function4(double x, double y)
{
    return x  * x + y * y * y - 1 - 8;
}

double func3dx(double x, double y) 
{
    return 3 * x * x;
}

double func3dy(double x, double y)
{
    return 2 * y;
}
double func4dx(double x, double y)
{
    return 2*x;
}
double func4dy(double x, double y)
{
    return 3 * y * y;
}
void inverse_matrix(double a[2][2])
{
    double det, aa;
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    aa = a[0][0];
    a[0][0] = a[1][1] / det;
    a[1][1] = aa / det;
    aa = a[0][1];
    a[0][1] = -a[0][1] / det;
    a[1][0] = -a[1][0] / det;
}

void SystemNewton(double x, double y)
{
    int i = 1;
    double a[2][2], dx, dy, b[2], norm;
    do
    {
        a[0][0] = func1dx(x, y);
        a[0][1] = func1dy(x, y);
        a[1][0] = func2dx(x, y);
        a[1][1] = func2dy(x, y);
        
        inverse_matrix(a);
        
        dx = -a[0][0] * function1(x, y) + -a[0][1] * function2(x, y);
        dy = -a[1][0] * function1(x, y) + -a[1][1] * function2(x, y);
        x = x + dx;
        y = y + dy;
        b[0] = function1(x, y);
        b[1] = function2(x, y);
        norm = sqrt(b[0] * b[0] + b[1] * b[1]);
        i++;
    } while (norm >= eps);
    cout << "Newton: " << endl;
    cout << x << " " << y << endl;
    cout << "Residual for the first equation: " << endl;
    cout << function1(x, y) << endl; 
    cout << "Residual for the second equation: " << endl;
    cout << function2(x, y) << endl;
}

void SystemNewton2(double x, double y)
{
    int i = 1;
    double a[2][2], dx, dy, b[2], norm;
    do
    {
        a[0][0] = func3dx(x, y);
        a[0][1] = func3dy(x, y);
        a[1][0] = func4dx(x, y);
        a[1][1] = func4dy(x, y);

        inverse_matrix(a);

        dx = -a[0][0] * function3(x, y) + -a[0][1] * function4(x, y);
        dy = -a[1][0] * function3(x, y) + -a[1][1] * function4(x, y);
        x = x + dx;
        y = y + dy;
        b[0] = function3(x, y);
        b[1] = function4(x, y);
        norm = sqrt(b[0] * b[0] + b[1] * b[1]);
        i++;
    } while (norm >= eps);
    cout << "Newton 2: " << endl;
    cout << x << " " << y << endl;
    cout << "Residual for the first equation: " << endl;
    cout << function3(x, y) << endl;
    cout << "Residual for the second equation: " << endl;
    cout << function4(x, y) << endl;
}

void ModifiedSystemNewton(double x, double y)
{
    int i = 1;
    double a[2][2], dx, dy, b[2], norm;
    double x0 = x;
    double y0 = y;
    a[0][0] = func1dx(x0, y0);
    a[0][1] = func1dy(x0, y0);
    a[1][0] = func2dx(x0, y0);
    a[1][1] = func2dy(x0, y0);

    inverse_matrix(a);
    do
    {
        dx = -a[0][0] * function1(x, y) + -a[0][1] * function2(x, y);
        dy = -a[1][0] * function1(x, y) + -a[1][1] * function2(x, y);
        x = x + dx;
        y = y + dy;
        b[0] = function1(x, y);
        b[1] = function2(x, y);
        norm = sqrt(b[0] * b[0] + b[1] * b[1]);
        i++;
    } while (norm >= eps);
    cout << "i = " << i << endl;
    cout << "Modified Newton: " << endl;
    cout << x << " " << y << endl;
    cout << "Residual for the first equation: " << endl;
    cout << function1(x, y) << endl;
    cout << "Residual for the second equation: " << endl;
    cout << function2(x, y) << endl;

}
void ModifiedSystemNewton2(double x, double y)
{
    int i = 1;
    double a[2][2], dx, dy, b[2], norm;
    double x0 = x;
    double y0 = y;
    dx = 0;
    dy = 0;
    a[0][0] = func3dx(x0, y0);
    a[0][1] = func3dy(x0, y0);
    a[1][0] = func4dx(x0, y0);
    a[1][1] = func4dy(x0, y0);
    double k = 0;
    inverse_matrix(a);
    do
    {
        cout << "x= " << x << endl;
        cout << "y = " << y << endl;
        dx = (a[0][0] * function3(x, y) + a[0][1] * function4(x, y));
        dy = -a[1][0] * function3(x, y) + -a[1][1] * function4(x, y);
        cout << "dx= " <<dx << endl;
        cout << "dy = " << dy << endl;
        double dx_better = dx * 1000.0;
        int dx_cut = dx_better;
        double dx_revert = dx_cut / 1000.0;

        double dy_better = dy * 1000.0;
        int dy_cut = dy_better;
        double dy_revert = dy_cut / 1000.0;

        k = x - dx_revert;
        cout << "k = " << k << endl;
        x = x - dx_revert;
        y = y + dy_revert;

        cout << "x= " << x << endl;
        cout << "y = " << y << endl;
        b[0] = function3(x, y);
        b[1] = function4(x, y);
        norm = sqrt(b[0] * b[0] + b[1] * b[1]);
        i++;

    } while (i<15);
    cout << "i = " << i<<  endl;
    cout << "Modified Newton 2: " << endl;
    cout << x << " " << y << endl;
    cout << "Residual for the first equation: " << endl;
    cout << function3(x, y) << endl;
    cout << "Residual for the second equation: " << endl;
    cout << function4(x, y) << endl;

}

void Relaxation(double x, double y)
{
    double a[2], dx, dy, x0 = 0, y0 = 0;

    int i = 0;
    do
    {
        a[0] = function1(x, y);
        a[1] = function2(x, y);

        dx = 0.1 * (-a[0]);
        dy = 0.1 * (-a[1]);
        x0 = x;
        y0 = x;
        x = x + dx;
        y = y + dy;

        i++;

    } while (i<100);

    cout << "Relaxation: " << endl;
    cout << x << " " << y << endl;
    cout << "Residual for the first equation: " << endl;
    cout << function1(x, y) << endl;
    cout << "Residual for the second equation: " << endl;
    cout << function2(x, y) << endl;
}

void Relaxation2(double x, double y)
{
    double a[2], dx, dy, x0 = 0, y0 = 0;

    int i = 0;
    do
    {
        a[0] = function3(x, y);
        a[1] = function4(x, y);

        dx = 0.1 * (-a[0]);
        dy = 0.1 * (-a[1]);
        x0 = x;
        y0 = x;
        x = x + dx;
        y = y + dy;

        i++;

    } while (i < 100);

    cout << "Relaxation 2: " << endl;
    cout << x << " " << y << endl;
    cout << "Residual for the first equation: " << endl;
    cout << function3(x, y) << endl;
    cout << "Residual for the second equation: " << endl;
    cout << function4(x, y) << endl;
}
int main()
{
    double x, y;
    x = 1, y = 1;
    double x0 = 1;
    double y0 = 1;
    SystemNewton(x, y);
    cout << endl;
    ModifiedSystemNewton(x, y);
    cout << endl;
    Relaxation(x, y);
    cout << "///////////////////////////////////////";
    SystemNewton2(x, y);
    cout << endl;
    ModifiedSystemNewton2(x0, y0);
    cout << endl;
    Relaxation2(x, y);
 
}



