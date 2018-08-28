#pragma once
#include<vector>
using namespace std;
class Matrix
{  
public:
    int n;
    int m;
    vector<vector<double>> a;
    ~Matrix();
    Matrix();
    Matrix(vector<double>, int, int);
    Matrix(vector<vector<double>>);
};
class Vector
{
public:
    int n;
    vector<double> a;
    ~Vector();
    Vector();
    Vector(vector<double>);
};
int main();
Vector gaussSolve(Matrix, Vector);
Vector gaussSeidelSolve(Matrix, Vector, Vector, int);
Vector tridiagonalSolve(Matrix, Vector);
void print(Vector);
void print(Matrix);
Vector randVect(int, int);
Matrix randTri(int, int);
Matrix RM(Matrix, int, double);
Matrix O(int, int);
Matrix I(int);
Matrix RM(int, int, double);
Matrix RP(int, int, int);
Matrix RA(int, int, int, double);
double abs(Vector);
double scalarProduct(Vector, Vector);
Matrix matrixMultiply(Matrix, Matrix);
Vector matrixVectorMultiply(Matrix, Vector);
Matrix matrixAdd(Matrix, Matrix);
Vector vectorAdd(Vector, Vector);
Matrix matrixScale(Matrix, double);
Vector vectorScale(Vector, double);
Matrix operator *(Matrix, Matrix);
Vector operator *(Matrix, Vector);
Matrix operator *(double, Matrix);
Vector operator *(double, Vector);
Matrix operator +(Matrix, Matrix);
Vector operator +(Vector, Vector);
Matrix operator -(Matrix, Matrix);
Vector operator -(Vector, Vector);
double operator *(Vector, Vector);
