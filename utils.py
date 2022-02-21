""" Author: Alberto Hernández López """
""" Utilidades que el programa principal necesita """
import  sys         as so
import  math        as math
import  numpy       as np

from    numpy     import *
from    colors    import *


# Regresa el indice de las filas y la 
# cantidad de zeros que contiene
# return idx_and_zeros -> tupla con indice y num de zeros
def how_many_zeros_you_have (matrix) :
    idx_and_zeros = []
    num_zeros       = 0
    (n, m) = matrix.shape
    for i in range(n - 1) : 
        column = matrix[i]
        for j in range(m - 1) : 
            element = column[j]
            if (fabs(element) < 1.0e-12) : 
                num_zeros += 1
        idx_and_zeros.append((i, num_zeros))

    return idx_and_zeros


# Revisa si el valor es cero
def is_a_zero (number) :
    return fabs(number) < 1.0e-12


# Cambia las filas de la matriz
def change_row (current_row, new_row, matrix) :
    aux_row = np.matrix(matrix[current_row, :])
    matrix[current_row, :] = matrix[new_row, :]
    matrix[new_row] = aux_row
    return matrix


# Revisa si hay filas en ceros
def get_rows_zeros_and_non_zero_indices (matrix) :
    (n, m)              = matrix.shape
    zero_indices        = []
    non_zero_indices    = []
    
    for k in range (n - 1) :
        (index, _) = get_non_zero_element(matrix[k, :])
        if (index == -1) :
            zero_indices.append(k)
        else :
            non_zero_indices.append(k)
    return (zero_indices, non_zero_indices)


# Mueve las filas llenas de ceros
def move_rows_zeros (zero_indices, non_zero_indices, matrix_a, matrix_b) :
    aux_i = len(non_zero_indices) - 1

    for index in zero_indices :         
        aux_row_a = matrix_a[index, :]
        matrix_a[index, :] = matrix_a[non_zero_indices[aux_i], :]
        matrix_a[non_zero_indices[aux_i], :] = aux_row_a

        aux_row_b = matrix_b[index, :]
        matrix_b[index, :] = matrix_b[non_zero_indices[aux_i], :]
        matrix_b[non_zero_indices[aux_i], :] = aux_row_b
        aux_i -= 1

    return matrix_a, matrix_b


# Normaliza la fila de la matriz
def normalize_row (factor, row) :
    return row / factor


# Obtiene el elemento no zero de la fila
# return (index, element)
def get_non_zero_element (row) : 
    for index in range(len(row)) : 
        element = row[index]
        if (not is_a_zero(element)) : 
            return (index, element)

    return (-1, 0)


# Obtiene el rango de la matriz
def get_matrix_range (matrix) :
    range_matrix    = 0
    (n, m)          = matrix.shape
    for i in range(n) :
        for j in range(m) :
            if (fabs(matrix[i, j]) != 0) :
                range_matrix += 1
                break
                
    return range_matrix


# Realiza una eliminacion final, en caso de necesitarla, de las ultimas dos filas
def check_final_rows (matrix_a, matrix_b) : 
    (n, m) = matrix_a.shape

    (index_n, value_a_n)     = get_non_zero_element(matrix_a[n - 1, :])
    (index_n_1, value_a_n_1) = get_non_zero_element(matrix_a[n - 2, :])

    value_b_n   = matrix_b[n - 1]
    value_b_n_1 = matrix_b[n - 2]

    if (index_n == index_n_1 and value_a_n == value_a_n_1 and value_b_n == value_b_n_1) :
        matrix_a[n - 1] = matrix_a[n - 1] - matrix_a[n - 2] * value_a_n_1
        matrix_b[n - 1] = matrix_b[n - 1] - matrix_b[n - 2] * value_a_n_1

    return (matrix_a, matrix_b)


# Obtiene la columna que no tiene zeros.
# return 
#   (i, j) -> coordenadas del valor que va a pivotear
#   value  -> el valor para pivotear
#   column -> columna para pivotear
def get_most_left_column (matrix_a) : 
    (n, m) = matrix_a.shape
    for i in range(m - 1) :
        column = matrix_a[:, i]
        for j in range(n) :
            value = column[j]
            if (not is_a_zero(value)) : 
                return ((j, i), value, column)
    
    return (None, None)


# Imprime las soluciones
def print_solutions (solutions, numb_excercise) : 
    print("Las soluciones del ejercicio", numb_excercise, "son: \n", solutions)


# Trunca un flotante a n digitos decimales
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper


# Imprime las ecuaciones que tienen soluciones infinitas
def print_infinite_solutions (matrix_a, matrix_b) :
    print("\n")
    (n, m) = matrix_a.shape
    solutions   = []
    solution    = ""
    for i in range(n) :
        for j in range(m) :
            if (not is_a_zero(matrix_a[i, j])) :
                value = matrix_a[i, j]
                if (value == 1.0) : 
                    solution += "X" + str(j)

                else :
                    solution += "(" + str(truncate(value, 2)) + ")" + "*X" + str(j)

                if j + 1 != m :
                    solution += " + "


        if (not solution == "") :
            solution += " = " + str(matrix_b[i]).lstrip('[').rstrip(']')
        print(f"{bcolors.OKGREEN}", solution)

        solution = ""        
    print("\n")

