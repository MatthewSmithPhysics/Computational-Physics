#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<ctime>
#include"linearsystems.h"

using namespace std;

//Matrix class methods
//---------------------------------------------------------
Matrix::~Matrix()
{
    this->n = 0;
    this->m = 0;
    this->a = {};
}
Matrix::Matrix()
{
    this->n = 0;
    this->m = 0;
    this->a = {};
}
Matrix::Matrix(vector<vector<double>> a)
{
    this->n = a.size();
    this->m = a[0].size();
    this->a = a;
}
Matrix::Matrix(vector<double> a, int n, int m)
{
    this->n = n;
    this->m = m;
    this->a = {};
    for(int i = 0; i < n; i++)
    {
        vector<double> row(m);
        for(int j = 0; j < m; j++) 
        {
            row[j] = a[m*i + j]; 
        }
        this->a.push_back(row);
    }
}
//-----------------------------------------------------------------
//Vector class methods
//-----------------------------------------------------------------
Vector::~Vector()
{
    this->n = 0;
    this->a = {};
}
Vector::Vector()
{
    this->n = 0;
    this->a = {};
}
Vector::Vector(vector<double> a)
{
    this->a = a;
    this->n = a.size();
}
//-----------------------------------------------------------------

void print(Matrix A)
{
    for(int i = 0; i < A.n; i++)
    {
        for(int j = 0; j < A.m; j++) cout << left << setw(20) << setprecision(16) << A.a[i][j];
        cout << endl;
    }
}

void print(Vector v)
{
    for(int i = 0; i < v.n; i++)
    {
        cout << v.a[i] << endl;
    }
}

int main()
{
    ofstream gaussOut;
    gaussOut.open("gaussout.txt");
    ofstream triOut;
    triOut.open("triout.txt");
    int M = 10000;
    gaussOut << "M = " << M << endl;
    gaussOut << left << setw(40) << "n" << setw(40) << "T(n)" << endl;
    triOut << "M = " << 1000*M  << endl;
    triOut << left << setw(40) << "n" << setw(40) << "T(n)" << endl;
    for(int n = 3; n < 10; n++)
    {
        cout << n << endl;
        long int gauss1 = time(NULL);
        for(int i = 0; i < M; i++) 
        {
            Matrix A = randTri(n, 100);
            Vector y = randVect(n, 100);
            Vector x = gaussSolve(A, y);
        }
        long int gauss2 = time(NULL);
        long int Tgauss = gauss2 - gauss1;
        gaussOut << left << setw(40) << n << setw(40) << Tgauss << endl;
        long int tri1 = time(NULL);
        for(int i = 0; i < 1000*M; i++)
        {
            Matrix A = randTri(n, 100);
            Vector y = randVect(n, 100);
            Vector x = tridiagonalSolve(A, y);
        }
        long int tri2 = time(NULL);
        long int Ttri = tri2 - tri1;
        triOut << left << setw(40) << n << setw(40) << Ttri << endl;
    }
}

/*
Performs Gaussian Elimination with Interchange to solve system of form Ax = b. Returns  vector x that solves the system. A should be n x n. 
*/
Vector gaussSolve(Matrix A, Vector b)
{
    int n = A.n;
    for(int j = 0; j < n; j++)
    {
        double max = 0.0;
        int k = j;
        for(int i = j; i < n; i++)
        {
            if(A.a[i][j]*A.a[i][j] >= max)
            {
                max = A.a[i][j]*A.a[i][j];
                k = i;
            }
        }
        A = RP(n, j, k)*A;
        b = RP(n, j, k)*b;
        double p = 1.0/A.a[j][j];
        A = RM(n, j, p)*A;
        b = RM(n, j, p)*b;
        for(int i = 0; i < n; i++)
        {
            if(i != j)
            {
                p = -A.a[i][j];
                A = RA(n, j, i, p)*A;
                b = RA(n, j, i, p)*b;
            }
        }
    }
    return b;
}

/*
Solves linear system Ax = b, by performing Gauss-Seidel Iteration for N steps. Returns vector x(N) that approximately satisfies the system. 
*/
Vector gaussSeidelSolve(Matrix A, Vector y, Vector x0, int N)
{
    int n = A.n;
    Vector x = x0;
    Matrix B = A;
    Vector z = y;
    for(int i = 0; i < n; i++)
    {
        double p = 1.0/B.a[i][i];
        B = RM(n, i, p)*B;
        z = RM(n, i, p)*z;
    }
    Matrix C = I(n) - B;
    
    for(int i = 0; i < N; i++)
    {
        x = C*x + z;
    }
    return x;
}

/*
Solves a linear system Ax = y, where A is a tridiagonal matrix. Returns vector x which solves the system. 
*/
Vector tridiagonalSolve(Matrix A, Vector y)
{
    int N = A.n;
    vector<double> a(N);
    vector<double> b(N);
    vector<double> c(N);
    for(int i = 0; i < N; i++)
    {
        b[i] = A.a[i][i];
        if(i == 0)
        {
            a[i] = 0.0;
            c[i] = A.a[i][i+1];
        }else if(i == N - 1)
        {
            a[i] = A.a[i][i-1];
            c[i] = 0.0;
        }else
        {
            a[i] = A.a[i][i-1];
            c[i] = A.a[i][i+1];
        }
    }
    vector<double> theta(N);
    vector<double> phi(N);
    theta[0] = -c[0]/b[0];
    phi[0] = y.a[0]/b[0];
    for(int i = 1; i < N; i++)
    {
        theta[i] = -c[i]/(a[i]*theta[i-1] + b[i]);
        phi[i] = (y.a[i] - a[i]*phi[i-1])/(a[i]*theta[i-1] + b[i]);
    }
    vector<double> x(N);
    x[N-1] = phi[N-1];
    for(int i = N - 2; i >= 0; i += -1) x[i] = theta[i]*x[i+1] + phi[i];
    Vector xf = Vector(x);
    return xf;
}

/*
Returns a randon n x n tridiagonal matrix, with integer elements between -M and M.
*/
Matrix randTri(int n, int M)
{
    Matrix A = O(n, n);
    A.a[0][0] = (rand() % 2*M) - M;
    A.a[0][1] = (rand() % 2*M) - M;
    A.a[n-1][n-2] = (rand() % 2*M) - M;
    A.a[n-1][n-1] = (rand() % 2*M) - M;
    for(int i = 1; i < n - 1; i++)
    {
        A.a[i][i-1] = (rand() % 2*M) - M;
        A.a[i][i] = (rand() % 2*M) - M;
        A.a[i][i+1] = (rand() % 2*M) - M;
    }
    return A;
}

/*
Returns a randon dimension n vector, with integer elements between -M and M
*/
Vector randVect(int n, int M)
{
    vector<double> x(n);
    for(int i = 0; i < n; i++) x[i] = (rand() % 2*M) - M;
    Vector v = Vector(x);
    return v;   
}

/*
Returns the n x m zero matrix
*/
Matrix O(int n, int m)
{
    vector<double> a(n*m);
    for(int i = 0; i < n*m; i++) a[i] = 0.0;
    Matrix A = Matrix(a, n, m);
    return A;
}

/*
Returns the n x n identity matrix
*/
Matrix I(int n)
{
    Matrix A = O(n, n);
    for(int i = 0; i < n; i++) A.a[i][i] = 1.0;
    return A;
}

/*
Returns the n x n matrix that, when applied to another matrix (m x n), multplies all elements in row k by p.
*/
Matrix RM(int n, int k, double p)
{
    Matrix A = I(n);
    A.a[k][k] = p;
    return A;
}

/*
Returns the n x n matrix that, when applied to another matrix (m x n), permutes rows i and j
*/
Matrix RP(int n, int i, int j)
{
    Matrix A = I(n);
    vector<double> row = A.a[i];
    A.a[i] = A.a[j];
    A.a[j] = row;
    return A;
}

/*
Returns the n x n matrix that, when applied to another matrix (m x n), adds p multiples of row i onto row j
*/
Matrix RA(int n, int i, int j, double p)
{
    Matrix A = I(n);
    A.a[j][i] = p;
    return A;
}

/*
Performs a row multilication on matrix A, kth row by p. Returns transformed matrix.
*/
Matrix RM(Matrix A, int k, double p)
{
    for(int j = 0; j < A.m; j++)
    {
        A.a[k][j] = A.a[k][j]*p;
    }
    return A;
}

/*
Computes modulus of vector v and returns the result.
*/
double abs(Vector v)
{
    double y = v*v;
    y = pow(y, 0.5);
    return y;
}

/*
Multiplies two matrices A and B, then returns the result. 
*/
Matrix matrixMultiply(Matrix A, Matrix B)
{
    Matrix C = Matrix();
    C.n = A.n;
    C.m = B.m;
    for(int i = 0; i < A.n; i++)
    {
        vector<double> row(B.m);
        for(int j = 0; j < B.m; j++)
        {
            double p = 0.0;
            for(int k = 0; k < A.m; k++)
            {
                p = p + A.a[i][k]*B.a[k][j];
            }
            row[j] = p;
        }
        C.a.push_back(row);
    }
    return C;
}

/*
Multiplies a matrix A and a vector v, then returns the result. 
*/
Vector matrixVectorMultiply(Matrix A, Vector v)
{
    Vector w = Vector();
    w.n = A.n;
    for(int i = 0; i < v.n; i++)
    {
        double p = 0.0;
        for(int j = 0; j < A.m; j++)
        {
            p = p + A.a[i][j]*v.a[j];
        }
        w.a.push_back(p);
    }
    return w;
}

/*
Adds matrix A to matrix B and returns the result
*/
Matrix matrixAdd(Matrix A, Matrix B)
{
    Matrix C = A;
    for(int i = 0; i < A.n; i++)
    {
        for(int j = 0; j < A.m; j++)
        {
            C.a[i][j] = C.a[i][j] + B.a[i][j];
        }
    }
    return C;
}

/*
Adds vector u to vector v and returns the result
*/
Vector vectorAdd(Vector u, Vector v)
{
    Vector w = u;
    for(int i = 0; i < u.n; i++)
    {
        w.a[i] = w.a[i] + v.a[i];
    }
    return w;
}

/*
Performs calculation of scalar product of u and v. Returns the result.
*/
double scalarProduct(Vector u, Vector v)
{
    double y = 0.0;
    for(int i = 0; i < u.n; i++)
    {
        y = y + u.a[i]*v.a[i];
    }
    return y;
}

/*
Multiplies matrix A by a scalar p and returns the result. 
*/
Matrix matrixScale(Matrix A, double p)
{
    Matrix B = A;
    for(int i = 0; i < A.n; i++)
    {
        for(int j = 0; j < A.m; j++)
        {
            B.a[i][j] = B.a[i][j]*p;
        }
    }
    return B;
}

/*
Multiplies vector v by scalar p and returns the result
*/
Vector vectorScale(Vector v, double p)
{
    Vector w = v;
    for(int i = 0; i < v.n; i++)
    {
        w.a[i] = w.a[i]*p;
    }
    return w;
}

Matrix operator *(Matrix A, Matrix B)
{
    Matrix C = matrixMultiply(A, B);
    return C;
}

Vector operator *(Matrix A, Vector v)
{
    Vector w = matrixVectorMultiply(A, v);
    return w;
}

Matrix operator +(Matrix A, Matrix B)
{
    Matrix C = matrixAdd(A, B);
    return C;
}

Vector operator +(Vector u, Vector v)
{
    Vector w = vectorAdd(u, v);
    return w;
}

Matrix operator *(double p, Matrix A)
{
    Matrix B = matrixScale(A, p);
    return B;
}

Vector operator *(double p, Vector v)
{
    Vector w = vectorScale(v, p);
    return w;
}

Matrix operator -(Matrix A, Matrix B)
{
    Matrix C = A + ((-1.0)*B);
    return C;
}

Vector operator -(Vector u, Vector v)
{
    Vector w = u + ((-1.0)*v);
    return w;
}

double operator *(Vector u, Vector v)
{
    double y = scalarProduct(u, v);
    return y;
}
