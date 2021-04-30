//cubic spline 

#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include <fstream>
#include "gnuplot.h"
using namespace std;

using vec = vector<double>;

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

    gp.sendLine("load 'f_values2.txt'");
  
    /*gp.sendLine("f(x) = (x > 0 && x <= 1) ? 0.2*x + 0.8*x**3 : 1/0");
    gp.sendLine("plot f(x)");
    gp.sendLine("h(x) = (x >= 1 && x <= 2) ? 2*x**3 - 3.6*x**2 + 3.8*x - 1.2 : 1/0");
    gp.sendLine("plot h(x)");
    gp.sendLine("g(x) = (x >= 2 && x < 4) ? -2.8*x**3 + 25.2*x**2 - 53.8*x + 37.2 : 1/0 ");
    gp.sendLine("plot g(x)"); 
    gp.sendLine("plot 'points_2.txt' w p ls 1 title 'data',  f(x), g(x), h(x)");*/

}
struct SplineSet {
    double a;
    double b;
    double c;
    double d;
    double x;
};

vector<SplineSet> spline(vec& x, vec& y)
{
    int n = x.size() - 1;
    vec a;
    a.insert(a.begin(), y.begin(), y.end());
    vec b(n);
    vec d(n);
    vec h;

    for (int i = 0; i < n; ++i)
        h.push_back(x[i + 1] - x[i]);

    vec alpha;
    alpha.push_back(0);
    for (int i = 1; i < n; ++i)
        alpha.push_back(3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1]);

    vec c(n + 1);
    vec l(n + 1);
    vec mu(n + 1);
    vec z(n + 1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n; ++i)
    {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for (int j = n - 1; j >= 0; --j)
    {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / 3 / h[j];
    }

    vector<SplineSet> output_set(n);
    for (int i = 0; i < n; ++i)
    {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
    }
    return output_set;
}

vector<double> standardPolynomial(const vector<double>& coef2, const vector<double>& x, int j)
{

    int N = x.size();
    vector<double> a(N, 0.0);                   // a[j] holds coefficient of x^j in final sum
    vector<double> p(N), prev(N);
    p[0] = 1;
    a[0] = coef2[0] * p[0];

    for (int i = 1; i < N; i++)
    {
        prev = p;
        if (x[j] != 0)
        {
            p[0] = -x[j] * prev[0];
        }
        else {
            p[0] = 0;
        }
        a[0] += coef2[i] * p[0];
        for (int k = 1; k <= i; k++)
        {
            p[k] = prev[k - 1] - x[j] * prev[k];
            a[k] += coef2[i] * p[k];
        }
    }
    return a;
}

int main()
{
    int n = 0;
    cout << "Enter number of points: ";
    cin >> n;
    vec x(n);
    vec y(n);
    for (int i = 0; i < x.size(); ++i)
    {
        x[i] = i;
        //y[i] = x[i]*x[i]*x[i];
        y[i] = (3.0*cos(2.0*x[i]));
    }

    std::ofstream out("points_2.txt");
    if (out.is_open())
    {
        for (int i = 0; i < x.size(); i++)
        {
            out << x[i] << " " << y[i] << std::endl;
        }
    }
    out.close();
    vector<SplineSet> cs = spline(x, y);
    vector<double> coef(n);
    vector<double> coef2(n);
    vector<vector<double>> allcoef;
    for (int j = 0; j < x.size() - 1; j++)
    {

        coef2[0] = cs[j].a;
        coef2[1] = cs[j].b;
        coef2[2] = cs[j].c;
        coef2[3] = cs[j].d;

        coef = standardPolynomial(coef2, x, j);
        allcoef.push_back(coef);
    }
    for (int i = 0; i < allcoef.size(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (j == 3)
            {
                cout << "(" << allcoef[i][j] << ")x**" << j << endl;
            }
            else
            {
                cout << "(" << allcoef[i][j] << ")x**" << j << "+";
            }
        }
        cout << endl;
    }
    std::ofstream out2("f_values2.txt");
    if (out2.is_open())
    {
        out2 << "plot 'points_2.txt' w p ls 1 title 'data', x**3 title 'main function',";
        for (int j = 0; j < allcoef.size(); j++)
        { 
            for (int i = 0; i < 4; i++)
            {
                if (i == 3 && i != 0)
                {
                    out2 << allcoef[j][i] << "*x**" << i;
                }
                else if (i == 0 && i != 3)
                {
                    out2 << allcoef[j][i] << "+";
                }
                else if (i == 0 && i == 3) {
                    out2 << allcoef[j][i];
                }
                else {
                    out2 << allcoef[j][i] << "*x**" << i << "+";
                }
            }
           // out2 << " x > " << x[j] << " &&" << " x <= " << x[j + 1] << " ";

            if(j!=allcoef.size() -1)
            {
                //out2 <<  " ,";
                out2 << ", ";
            }
            else {
                out2 << "";
            }

        }
    }
    out2.close();


    for (int i = 0; i < x.size() - 1; i++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            std::cout << allcoef[i][j] << "\t";
        }
        std::cout << "\n";
    }



    drawGnuplot();
}