#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include <fstream>
#include <string>
#include <cstring>
#include "gnuplot.h"


using namespace std;

double func(double x)
{
    return sin(x);
}
void drawGnuplot()
{
    GnuplotPipe gp;
    gp.sendLine("set border linewidth 1.5");
    gp.sendLine("set style line 1 lc rgb '#0060ad' pt 7 ps 1.5 lt 1 lw 2 # -- - blue");
    gp.sendLine("unset key");
    gp.sendLine("set style line 11 lc rgb '#808080' lt 1");
    gp.sendLine("set border 3 back ls 1");
    gp.sendLine("set border 3 back ls 1");
    gp.sendLine("set border 3 back ls 1");
    gp.sendLine("set tics nomirror out scale 0.75");
    gp.sendLine("set style line 12 lc rgb'#808080' lt 0 lw 1");
    gp.sendLine("set grid back ls 12");
    gp.sendLine("set xrange[-10:20]");
    gp.sendLine("set yrange[-10:30]");
    gp.sendLine("set key invert");
    gp.sendLine("set key right top");
   
    gp.sendLine("load 'f_values.txt'");
    
}


void Gauss(vector<double> points, vector<double> f_values)
{
    int N = points.size();
    vector<vector<double>> matrix(N, std::vector<double>(N + 1, 0));
    //double matrix[N][N + 1];
    for (int m1 = 0; m1 < N; m1++)
    {
        for (int m2 = 0; m2 < N; m2++)
        {
            matrix[m1][m2] = pow(points[m1], m2);
        }
    }
    for (int m1 = 0; m1 < N; m1++)
    {
        for (int m2 = N; m2 < N + 1; m2++)
        {
            matrix[m1][m2] = f_values[m1];
        }
    }
    cout << "\n ---------------------------------\n";
    cout << "\n  Matrix is:\n";
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
            cout << setw(8) << setprecision(4) << matrix[i][j];
        cout << endl;
    }
    double temp, s;
    //double x[N];
    vector<double> x(N);
    //make above matrix upper triangular Matrix

    for (int j = 0; j < N - 1; j++)
    {
        for (int i = j + 1; i < N; i++)
        {
            temp = matrix[i][j] / matrix[j][j];

            for (int k = 0; k < N + 1; k++)
                matrix[i][k] -= matrix[j][k] * temp;
        }
    }


    cout << "\n ---------------------------------\n";
    cout << "\n Upper Triangular Matrix is:\n";
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
            cout << setw(8) << setprecision(4) << matrix[i][j];
        cout << endl;
    }


    cout << "\n ---------------------------------\n";

    for (int i = N - 1; i >= 0; i--)
    {
        s = 0;
        for (int j = i + 1; j < N; j++)
            s += matrix[i][j] * x[j];
        x[i] = (matrix[i][N] - s) / matrix[i][i];
    }


    /*cout << "\n Solution is:\n";
    for (int i = 0; i < N; i++)
        cout << "x[" << setw(3) << i + 1 << "]=" << setw(7) << setprecision(4) << x[i] << endl;*/
    cout << "Gaussian Polinomial is: " << endl;
    std::ofstream out("f_values.txt");
    if (out.is_open())
    {
        out << "plot 'points.txt' w p ls 1 title 'data', ";
        for (int i = 0; i < N; i++)
        {
            if (i == N - 1 && i!=0)
            {
                out  << x[i] << "*x**" << i ;
            }
            else if (i == 0 && i != N - 1)
            {

                out << x[i] << "+";
            }
            else if (i == 0 && i == N - 1) {
                out << x[i];
            }
            
            else {

                out  << x[i] << "*x**" << i << "+";
            }
        }
        out << " title 'Gauss'";
    }
    out.close();
   
    for (int i = 0; i < N; i++)
    {
        if (i == N - 1)
        {
            cout << "(" << x[i] << ")x**" << i << endl;
        }
        else {
            cout << "(" << x[i] << ")x**" << i << "+";
        }
    }

}
void Lagrange(vector<double> points, vector<double> f_values)
{

    auto res = new vector<double>(points.size(), 0); //final coefficients
    for (int i = 0; i < points.size(); i++)
    {
        vector<double> tmpcoeffs(points.size(), 0);

        // Start with a constant polynomial
        tmpcoeffs[0] = f_values[i];
        double prod = 1;
        for (auto point : points)
        {
            if (points[i] == point) continue;
            prod *= points[i] - point;
            double precedent = 0;
            for (auto resptr = tmpcoeffs.begin(); resptr < tmpcoeffs.end(); resptr++)
            {
                // Compute the new coefficient of X^i based on
                // the old coefficients of X^(i-1) and X^i
                double newres = (*resptr) * (-point) + precedent;
                precedent = *resptr;
                *resptr = newres;
            }


        }
        /*transform(res->begin(), res->end(),
            tmpcoeffs.begin(),
            res->begin(),
            [=](double oldcoeff, double add) {return oldcoeff + add / prod; }
        );*/
        for (int i = 0; i < points.size(); i++)
        {
            double num = tmpcoeffs[i] / prod;
            double value = res->operator[](i);
            res->at(i) = num + value;
        }

    }

    cout << "Lagrange Polinomial is: " << endl;
    std::ofstream out("f_values.txt", std::ios::app);
    if (out.is_open())
    {
        out << ", ";
        for (int i = 0; i < points.size(); i++)
        {
            if (i == points.size() - 1 && i != 0)
            {
                out << res->at(i) << "*x**" << i ;
            }
            else if (i == 0 && i != points.size() - 1)
            {

                out << res->at(i) << "+";
            }
            else if (i == 0 && i == points.size() - 1) {
                out << res->at(i) << endl;
            }

            else {

                out << res->at(i) << "*x**" << i << "+";
            }
        }
        out << " title 'Lagrange'";
    }
    out.close();

    for (int i = 0; i < points.size(); i++)
    {
        if (i == points.size() - 1)
        {
            cout << "(" << res->at(i) << ")x**" << i << endl;
        }
        else {
            cout << "(" << res->at(i) << ")x**" << i << "+";
        }
    }



}
void Newton(const vector<double>& x, vector<double> f)          
{
    int N = x.size();
    vector<double> c(N), temp(N);

    c[0] = f[0];
    for (int i = 1; i < N; i++)       // Compute ith divided differences
    {
        for (int j = 0; j < N - i; j++) temp[j] = (f[j + 1] - f[j]) / (x[j + i] - x[j]);
        f = temp;
        c[i] = f[0];
    }
    

    
    vector<double> a(N, 0.0);                   // a[j] holds coefficient of x^j in final sum
    vector<double> p(N), prev(N);                 // At the ith step, p[ ] is the ith polynomial of Newton factors

    p[0] = 1;
    a[0] = c[0] * p[0];
    for (int i = 1; i < N; i++)
    {
        prev = p;
        p[0] = -x[i - 1] * prev[0];
        a[0] += c[i] * p[0];
        for (int j = 1; j <= i; j++)
        {
            p[j] = prev[j - 1] - x[i - 1] * prev[j];
            a[j] += c[i] * p[j];
        }
    }

    cout << "Newton Polinomial is: " << endl;
    std::ofstream out("f_values.txt", std::ios::app);
    if (out.is_open())
    {
        out << ", ";
        for (int i = 0; i < x.size(); i++)
        {
            if (i == x.size() - 1 && i != 0)
            {
                out << a[i] << "*x**" << i;
            }
            else if (i == 0 && i != x.size() - 1)
            {

                out << a[i] << "+";
            }
            else if (i == 0 && i == a.size() - 1) {
                out << a[i] << endl;
            }

            else {

                out << a[i] << "*x**" << i << "+";
            }
        }
        out << "title 'Newton'";
    }
    out.close();

    for (int i = 0; i < x.size(); i++)
    {
        if (i == x.size() - 1)
        {
            cout << "(" << a[i] << ")x**" << i << endl;
        }
        else {
            cout << "(" << a[i] << ")x**" << i << "+";
        }
    }

   
}




int main()
{
    double a, b;
    cout << "Enter a:\n";
    cin >> a;
    cout << "Enter b:\n";
    cin >> b;
    int num_points;
    cout << "Enter the number of points:\n";
    cin >> num_points;
    vector<double> points;
    double interval = (b - a) / num_points;
    for (int i = 1; i <= num_points; i++)
    {
        double point = (b - a) * i / num_points;
        points.push_back(point);
    }
 
    vector<double> f_values;
    for (int i = 0; i < num_points; i++)
    {
        f_values.push_back(func(points[i]));
    }
    std::ofstream out("points.txt");
    if (out.is_open())
    {
        for (int i = 0; i < points.size(); i++)
        {
            out << points[i] << " " << f_values[i] << std::endl;
        }
    }
    out.close();
    Gauss(points, f_values);
    Lagrange(points, f_values);
    Newton(points, f_values);


    drawGnuplot();

}
