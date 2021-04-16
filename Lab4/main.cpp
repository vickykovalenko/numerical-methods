#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include <fstream>

using namespace std;

float func(float x)
{
    return sin(x);
}
void Gauss(vector<float> points, vector<float> f_values, int N)
{
    vector<vector<float>> matrix(N, std::vector<float>(N + 1, 0));
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
    float temp, s;
    //double x[N];
    vector<float> x(N);
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


    //print values of x,y,z

    cout << "\n The Solution is:\n";
    for (int i = 0; i < N; i++)
        cout << "x[" << setw(3) << i + 1 << "]=" << setw(7) << setprecision(4) << x[i] << endl;
    cout << "Polinomial is: " << endl;
    std::ofstream out("D:\\f_values.txt");
    if (out.is_open())
    {
        for (int i = 0; i < N; i++)
        {
            if (i == N - 1 && i!=0)
            {
                out  << x[i] << "*x**" << i << endl;
            }
            else if (i == 0 && i != N - 1)
            {

                out << x[i] << "+";
            }
            else if (i == 0 && i == N - 1) {
                out << x[i] << endl;
            }
            
            else {

                out  << x[i] << "*x**" << i << "+";
            }
        }
    }
    out.close();
    for (int i = 0; i < N; i++)
    {
        if (i == N - 1)
        {
            cout << "(" << x[i] << ")x^" << i << endl;
        }
        else {
            cout << "(" << x[i] << ")x^" << i << "+";
        }
    }



}

int main()
{
    float a, b;
    cout << "Enter a:\n";
    cin >> a;
    cout << "Enter b:\n";
    cin >> b;
    int num_points;
    cout << "Enter number of points:\n";
    cin >> num_points;
    vector<float> points;
    float interval = (b - a) / num_points;
    while (a < b)
    {
        float point = a + interval;
        points.push_back(point);
        a += interval;
    }
    vector<float> f_values;
    for (int i = 0; i < num_points; i++)
    {
        f_values.push_back(func(points[i]));
    }
    std::ofstream out("D:\\operations.txt");
    if (out.is_open())
    {
        for (int i = 0; i < points.size(); i++)
        {
            out << points[i] << " " << f_values[i] << std::endl;
        }
    }
    out.close();
    Gauss(points, f_values, num_points);

}
