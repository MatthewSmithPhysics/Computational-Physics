import numpy
from numpy import transpose, cos, sin, log, sqrt

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

def sym2tri(A):
    B = numpy.copy(A)
    n = len(A)
    M = I(n)
    for i in range(n, 2, -1):
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

def herm(A):
    B = numpy.copy(A)
    B = transpose(B).conj()
    return B

def cnorm(A):
    y = (herm(A)@A)[0, 0]
    y = y**0.5
    return y

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