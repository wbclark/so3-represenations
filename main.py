import numpy as np
from cmath import sqrt

def spin_to_dim(spin):
    return int(1+2*spin)

def kronecker_delta(a, b):
    if a == b:
        return 1.0
    else:
        return 0.0

def Jx(spin):
    dim = spin_to_dim(spin)
    matrix_real = np.zeros( (dim, dim) )
    matrix_imag = np.zeros( (dim, dim) )
    for a in range(dim):
        for b in range(dim):
            row_index, column_index = a+1, b+1
            scale_factor = 0.5*sqrt((spin+1)*(column_index+row_index-1)-column_index*row_index)
            presence_factor = kronecker_delta(row_index, column_index+1) + kronecker_delta(row_index+1, column_index)
            matrix_real[a][b] = scale_factor.real*presence_factor
            matrix_imag[a][b] = scale_factor.imag*presence_factor
    return matrix_real+(0+1j)*matrix_imag

def Jy(spin):
    dim = spin_to_dim(spin)
    matrix_real = np.zeros( (dim, dim) )
    matrix_imag = np.zeros( (dim, dim) )
    for a in range(dim):
        for b in range(dim):
            row_index, column_index = a+1, b+1
            scale_factor = (1/(0+2j))*sqrt((spin+1)*(column_index+row_index-1)-column_index*row_index)
            presence_factor = kronecker_delta(row_index, column_index+1) - kronecker_delta(row_index+1, column_index)
            matrix_real[a][b] = scale_factor.real*presence_factor
            matrix_imag[a][b] = scale_factor.imag*presence_factor
    return matrix_real+(0+1j)*matrix_imag

def Jz(spin):
    dim = spin_to_dim(spin)
    matrix = np.zeros( (dim, dim) )
    for a in range(dim):
        for b in range(dim):
            row_index, column_index = a+1, b+1
            scale_factor = spin+1-column_index
            matrix[a][b] = scale_factor*kronecker_delta(row_index, column_index)
    return matrix

def print_generators(spin):
    print("printing generators for spin: %s" % spin)
    print()
    print("Jx(%s)" % spin)
    print()
    print(Jx(spin))
    print()
    print("Jy(%s)" % spin)
    print()
    print(Jy(spin))
    print()
    print("Jz(%s)" % spin)
    print()
    print(Jz(spin))

print_generators(1.5)
