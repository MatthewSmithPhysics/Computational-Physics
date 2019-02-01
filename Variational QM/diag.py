

import numpy
import numpy.random as random
import matplotlib.pyplot as pyplot
import time
from numpy import transpose, cos, sin, log, sqrt
import sys


def I(n):
    M = numpy.zeros((n, n), dtype = float)
    for i in range(n):
        M[i, i] = 1.0
    return M

def basis(n, k):
    E = numpy.zeros((n, 1), dtype = float)
    E[k, 0] = 1.0
    return E

def norm(M):
    N = numpy.copy(M)
    y = (transpose(N)@N)[0, 0]
    y = y**0.5
    return y

def isTri(A):
    res = True
    n = len(A)
    for j in range(2, n):
        if(A[0, j] != 0.0): return False
    for j in range(n - 2):
        if(A[n - 1, j] != 0.0): return False
    for i in range(1, n - 1):
        for j in range(i - 1):
            if(A[i, j] != 0.0): return False
        for j in range(i + 2, n):
            if(A[i, j] != 0.0): return False
    return True

def sym2tri(A):
    B = numpy.copy(A)
    n = len(A)
    M = I(n)
    for i in range(n, 2, -1):
        if(isTri(B)): return B, M
        C = numpy.zeros((n, 1), dtype = float)
        C[0:(i - 1), 0] = B[0:(i - 1), i - 1]
        U = C - norm(C)*basis(n, i - 2)
        P = I(n) - 2.0*(U@transpose(U))/(norm(U)**2.0)
        B = P@(B@P)
        M = M@P
        
    return B, M


def gramm_schmidt(A):
    E = numpy.copy(A)
    n = len(A)
    for k in range(n):
        Ek = E[:, k]
        for i in range(0, k):
            Ek = Ek - numpy.dot(A[:, k], E[:, i])*E[:, i]
        Ek = numpy.dot(Ek, Ek)**(-0.5)*Ek
        E[:, k] = Ek
    return E

def QR(A):
    Q = gramm_schmidt(A)
    n = len(A)
    R = numpy.zeros((n, n), dtype = float)
    for i in range(n):
        for j in range(i, n):
            R[i, j] = numpy.dot(A[:, j], Q[:, i])
    return Q, R




def dense_sym_diagonalise_QR(A, N):
    n = len(A)
    B = numpy.copy(A)
    
    B, M = sym2tri(B)
    
    for i in range(N):
        Q, R = QR(B)
        B = R@Q
        M = M@Q
        
    D = numpy.copy(B)
    for i in range(n):
        for j in range(n):
            if(not i == j): D[i, j] = 0.0
                
    return D, M




def phi(A, k, t):
    y = 0.0
    if(k == 0): y = 1.0
    elif(k == 1): y = A[0, 0] - t
    else:
        y = (A[k - 1, k - 1] - t)*phi(A, k - 1, t) - A[k - 2, k - 1]**2.0*phi(A, k - 2, t)
    return y

def diff(a, b):
    if(a >= 0 and b <= 0): return True
    elif(a <= 0 and b >= 0): return True
    else: return False
    
def sigma(A, c):
    n = len(A)
    s = 0
    p = phi(A, 0, c)
    for i in range(n):
        q = phi(A, i + 1, c)
        if(diff(p, q)): s = s + 1
        p = q
    return s

def bisect_solve(A, a, b, M):
    n = len(A)
    p = phi(A, n, a)
    q = phi(A, n, b)
    for i in range(M):
        c = 0.5*(a + b)
        r = phi(A, n, c)
        if(diff(p, r)): 
            b = c
            q = r
        else: 
            a = c
            p = r
    x = 0.5*(a + b)
    return x

def sturm_solve(A, a, b, M):
    n = len(A)
    p = sigma(A, a)
    q = sigma(A, b)
    N = abs(p - q)
    sol = []
    if(N == 0):
        return []
    elif(N == 1):
        x = bisect_solve(A, a, b, M)
        sol.append(x)
    else:
        c = 0.5*(a + b)
        X = sturm_solve(A, a, c, M)
        Y = sturm_solve(A, c, b, M)
        sol = sol + X + Y
    return sol


def det(A):
    y = 0.0
    n = len(A)
    if(n == 1): return A[0, 0]
    for i in range(n):
        B = numpy.zeros((n - 1, n - 1), dtype = complex)
        for i in range(0, i):
            B[1:, i] = A[1:, i]
        for i in range(i + 1, n):
            B[1:, i - 1] = A[1:, i]
        y = y + A[0, i]*det(B)
    return y

def gauss_inverse(A):
    n = len(A)
    B = numpy.copy(A)
    C = I(n)
    for i in range(n):
        if(B[i, i] == 0.0):
            for j in range(i + 1, n):
                if(B[j, i] != 0.0): 
                    row = B[i, :]
                    B[i, :] = B[j, :]
                    B[j, :] = row
                    row = C[i, :]
                    C[i, :] = C[j, :]
                    C[j, :] = row
    
        m = B[i, i]
        B[i, :] = B[i, :]/m
        C[i, :] = C[i, :]/m
        for j in range(n):
            if(j != i):
                k = B[j, i]
                B[j, :] = B[j, :] - k*B[i, :]
                C[j, :] = C[j, :] - k*C[i, :]
    return C

def dense_sym_diagonalise_char(A, a, b, M):
    n = len(A)
    lam = sturm_solve(A, a, b, M)
    P = numpy.zeros((n, n), dtype = complex)
    for i in range(n):
        y = lam[i]
        ks = []


def herm(A):
    B = numpy.copy(A)
    B = transpose(B).conj()
    return B

def dense_herm_diagonalise_QR(H, N):
    n = len(H)
    A = numpy.real(H)
    B = numpy.imag(H)
    C = numpy.zeros((2*n, 2*n), dtype = float)
    C[0:n, 0:n] = A[:, :]
    C[0:n, n:] = -B[:, :]
    C[n:, 0:n] = B[:, :]
    C[n:, n:] = A[:, :]
    Dc, Mc = dense_sym_diagonalise_QR(C, N)
    D = numpy.zeros((n, n), dtype = complex)
    M = numpy.zeros((n, n), dtype = complex)

    TRACKER = numpy.copy(Dc) # matrix to aid selection of eigenvalues
    for j in range(n):
        pri = []
        for i in range(2*n): pri.append(abs(TRACKER[i, i]))
        k = pri.index(max(pri))
        D[j, j] = Dc[k, k]
        for i in range(n): M[i, j] = Mc[i, k] + 1.0j*Mc[n + i, k]
        TRACKER = TRACKER - Dc[k, k]*I(2*n)

    return D, M


