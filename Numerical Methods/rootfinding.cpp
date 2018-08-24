#include<iostream>
#include<cmath>
#include<vector>
#include<iomanip>
#include"rootfinding.h"

using namespace std;

/*
Equation to solve: x**5 - 2*x**2 - 3 = 0
Solution: c
*/

/*
For equation in form f(x) = 0
*/
double f(double x)
{
    double y = dpow(x, 5.0) - 2.0*dpow(x, 2.0) - 3.0;
    return y;
}


/*
Derivative of f
*/
double df(double x)
{
    double y = 5.0*dpow(x, 4.0) - 4.0*x;
    return y;
}

/*
For equation in form x = g(x)
*/
double g(double x)
{
    double y = dpow(2.0*dpow(x, 2.0) + 3.0, 0.2);
    return y;
}

int main()
{
    double a = newtonRaphson(1.5, 100, false, 0.0);
    cout << "Equation : x**5 - 2*x**2 - 3 = 0"<< endl;
    cout << left << setw(4) << "n" << setw(20) << "x(n)" << "|x(n) - a|" << endl;
    for(int i = 0; i < 100; i++)
    {
        vector<double> c = binaryChopping(-20.0, 10.0, i);
        cout << left << setw(4) << i << setw(20) << setprecision(16) << c[0] << abs(c[0] - a) << endl;
    }
    return 0;    
}

/*
Uses the rearranged form of the equation, x = g(x), to approach the solution to the equation. This is by recurrence relationship x(n + 1) = g[x(n)]. Function returns x(N). 
*/
double rearrangedIteration(double x0, int N)
{
    double x = x0;
    for(int i = 0; i < N; i++)
    {
        x = g(x);
    }
    return x;
}

/*
Uses the Newton-Raphson recurrence relationship x(n + 1) = x(n) - f[x(n)]/f'[x(n)] to approach the solution to the equation. numerical = true if numerical differentiation (with step size dx) is to be used. Function returns x(N)
*/
double newtonRaphson(double x0, int N, bool numerical, double dx)
{  
    double x = x0;
    for(int i = 0; i < N; i++)
    {
        double m;
        if(numerical) m = (f(x + dx) - f(x))/dx;
        else m = df(x);
        x = x - f(x)/m;
    }
    return x;
}

/*
At each step, the line segment connecting (A(n), f[A(n)]) and (B(n), f[B(n)]) is used to find a new estimate for the solution, by interpolation. This estimate then replaces one of the line's points for the next iteration step. Function returns a vector containing the final estimate and an error margin. 
*/
double interpolation(double A0, double B0, int N)
{
    double c;
    if(abs(f(B0)) > abs(f(B0) - f(A0)))
    {
        cout << "The interval [" << A0 << ", " << B0 << "] does not contain the solution." << endl;
        return c;
    }
    if(f(A0) == 0.0)
    {
        c = A0;
    }else if(f(B0) == 0.0)
    {
        c = B0;
    }else
    {
        double A = A0;
        double B = B0;
        double x = (A0*f(B0) - B0*f(A0))/(f(B0) - f(A0));
        int pos; //indicates which side of interval is +ve
        if(f(A0) > 0) pos = 0;
        else pos = 1;
        for(int i = 0; i < N; i++)
        {
            if(f(x) == 0.0)
            {
                c = x;
                return c;
            }
            if(f(x) > 0.0)
            {
                if(pos == 0) A = x;
                else B = x;
            }else
            {
                if(pos == 0) B = x;
                else A = x;
            }
            x = (A*f(B) - B*f(A))/(f(B) - f(A));
        }
        c = x;
    }
    return c;
}

/*
At each step, the midpoint (x, f(x)) of the line segment connecting (A(n), f[A(n)]) and (B(n), f[B(n)]) is chosen to replace A or B, depending on the sign of f(x). This has the effect of repeatedly halving the interval containing the solution to the equation. Function returns the final estimate and an error margin.
*/ 
vector<double> binaryChopping(double A0, double B0, int N)
{
    vector<double> c = {0.0, 0.0};
    if(abs(f(B0)) > abs(f(B0) - f(A0)))
    {
        cout << "The interval [" << A0 << ", " << B0 << "] does not contain the solution." << endl;
        return c;
    }
    if(f(A0) == 0.0)
    {
        c[0] = A0;
        c[1] = 0.0;
    }else if(f(B0) == 0.0)
    {
        c[0] = B0;
        c[1] = 0.0;
    }else
    {
        double A = A0;
        double B = B0;
        double x = 0.5*(A0 + B0);
        int pos; //indicates which of A or B is +ve
        if(f(A0) > 0.0) pos = 0;
        else pos = 1;
        for(int i = 0; i < N; i++)
        {
            if(f(x) == 0.0)
            {
                c[0] = x;
                c[1] = 0.0;
                return c;
            }
            if(f(x) > 0.0)
            {
                if(pos == 0) A = x;
                else B = x;
            }else
            {
                if(pos == 0) B = x;
                else A = x;
            }
            x = 0.5*(A + B);
        }
        c[0] = x;
        if(abs(x - A) > abs(B - x)) c[1] = abs(x - A);
        else c[1] = abs(B - x);
    }
    return c;
}
/*
Custom function to get the double result from cmath pow.
*/
double dpow(double x, double p)
{
    double y = pow(x, p);
    return y;
}
