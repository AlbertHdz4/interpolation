"""" Alberto Hernández López """
"""" Interpolacion de splines lineal 
    y cubico para funciones 1D """

import matplotlib.pyplot    as plt
import numpy                as np
import sympy                as sym
import bisect
import math
import time

from colors import bcolors


# Spline lineal para la interpolacion de funciones 1D
def linear_spline (x_i, f_i) :
    """
    Interpola funciones 1D usando un spline lineal
    x_i : Coordenadas x a ser evaluadas en f(x) y a ser interpolados 
    f_i : Valores o salidas de la funcion

    Se deben de sumar polinomio de la forma, 
    Pi(x) = ax + b
    """
    n = len(x_i) # Longitud de nuestro espacio
    x = sym.Symbol('x') # Indicamos que la variable x es simbolica y sera evaluada en los polinomios

    all_polys   = [] # Lista para ir almacenando los polinomios de cada subdominio
    position    = 1 # iterador para ir construyendo pendientes y polinomios

    # Construimos los polinomios 
    for position in range(1, n) :
        m = (f_i[position] - f_i[position - 1]) / (x_i[position] - x_i[position - 1]) # Calculamos la pendiente entre m = y[i] - y[i - 1] / x_i[i] - x_i[i - 1]
        poly_i = f_i[position - 1] 
        poly_i += m * (x - x_i[position - 1])
        all_polys.append(poly_i) # Agregamos el polinomio i-esimo para despues graficarlos

    """ Se imprimen todos los polinomios por segmento """
    # for i in range(n - 1) :
    #     print(all_polys[i])

    """ Subdominios y subrangos de cada polinomio """
    x_dom   = np.array([])
    y_range = np.array([])

    """ Guardamos los polinomios junto con su rango en el que estan definidos """
    for i in range(1, n) :
        lower_dom_i     = x_i[i - 1]
        upperr_dom_i    = x_i[i]
        sub_dom = np.linspace(lower_dom_i, upperr_dom_i, 2) # Muestras para graficar
        poly_i = all_polys[i - 1]
        poly_x_i = sym.lambdify(x, poly_i) # Expresion de polinomio a evaluar
        poly_y_i = poly_x_i(sub_dom) # Evaluando el polinomio
        x_dom = np.concatenate((x_dom, sub_dom))
        y_range = np.concatenate((y_range, poly_y_i))
        

    """ Graficamos spline lineal """
    plt.figure('Spline lineal')
    plt.title('Spline lineal')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.plot(x_i, f_i, 'o', color = 'red', label = 'Puntos de la muestra')
    plt.plot(x_dom, y_range, label = 'Spline lineal')
    plt.legend()
    plt.show()
    return poly_i
    

# Spline cubico para interpolacion de funciones de 1D
def get_cubic_spline (x_i, f_i):
    initial_time = time.time()

    n       = len(x_i) # Obtenemos la longitud de nuestros datos
    h_i     = get_h_diff(x_i) # obtenemos los cambios para obtener la 
    A, B, C = create_tridiagonal_matrix(n, h_i) # creamos la matriz triagonal
    D       = get_d_values(n, h_i, f_i) # obtenemos los dj valores para obtener las soluciones de la matriz tridiagonal
    M       = triadiagonal_matrix_solution(A, B, C, D) # obtenemos las soluciones de la matriz

    coefficients = [[(M[i + 1] - M[i]) * h_i[i] * h_i[i] / 6, M[i] * h_i[i] * h_i[i] / 2,\
                    (f_i[i + 1] - f_i[i] - (M[i + 1] + 2 * M[i]) * h_i[i] * h_i[i] / 6),\
                    f_i[i]] for i in range(n-1)]

    # print("Coefficients: ", coefficients)

    def spline(val):
        idx = min(bisect.bisect(x_i, val) - 1, n - 2)
        z = (val - x_i[idx]) / h_i[idx]
        C = coefficients[idx]
        return (((C[0] * z) + C[1]) * z + C[2]) * z + C[3]

    cubic_spline_y  = [spline(y) for y in x_i]
    
    """ Obtenemos los tiempos finales que tardo el algoritmo """
    final_time  = time.time()
    diff        = final_time - initial_time
    print('Spline cubico: El tiempo que tardo en solucionar el problema con ', n, ' puntos fue de ', diff, ' segundos')

    """ Graficamos spline cubico """
    plt.figure('Spline cubico')
    plt.title('Spline cubico')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.plot(x_i, f_i, 'o', color = 'red', label = 'Puntos de la muestra')
    plt.plot(x_i, cubic_spline_y, label = 'Spline cubico')
    plt.legend()
    plt.show()
    return spline


def create_tridiagonal_matrix (n, h) :
    """
    Interpola funciones 1D usando un spline cubico
    x_i     : Coordenadas x a ser evaluadas en f(x) y a ser interpolados 
    f_i     : Valores o salidas de la funcion
    spline  : Regresa los valores evaluados del spline cubico
    """
    A = np.zeros(n - 1)
    B = np.zeros(n)
    C = np.zeros(n - 1)

    for i in range(n - 2) : 
        A[i] = h[i] / (h[i] + h[i + 1])

    A[n - 2] = 0
    
    for j in range(n) : 
        B[j] = 2 
   
    C = [0] + [h[i + 1] / (h[i] + h[i + 1]) for i in range(n - 2)]
    return A, B, C


# Obtenemos las h[j] = x[j + 1] - x[j] de todos lo datos 
def get_h_diff (x) :
    h = np.zeros(len(x))
    for i in range(len(x) - 1) :
        h[i] = x[i + 1] - x[i]
    return h


# Obtenemos los d[j] valores
# d = d[0] + sum(6((y[j + 1] - y[j]/ h[i]) - (y[j] - y[i - 1])/h[j - 1]) / h[j] + h[j - 1]) + d[n]
# donde d[0] = d[n] = 0 
def get_d_values (n, h, y) :
    return [0] + [6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]) / (h[i] + h[i-1]) for i in range(1, n - 1)] + [0]


## Obtenemos las soluciones de la matriz triadiagonal
def triadiagonal_matrix_solution (A, B, C, D):
    c_p = C + [0]
    d_p = np.zeros(len(B), float)
    X   = np.zeros(len(B), float)

    c_p[0] = C[0] / B[0]
    d_p[0] = D[0] / B[0]

    for i in range(1, len(B)):
        c_p[i] = c_p[i] / (B[i] - c_p[i - 1] * A[i - 1])
        d_p[i] = (D[i] - d_p[i - 1] * A[i - 1]) / (B[i] - c_p[i - 1] * A[i - 1])

    # Realizamos la solucion hacia atras, de la ultima variable
    X[-1] = d_p[-1] # Primer elemento de la solucion

    # Vamos solucionando el sistema de adelante hacia atras
    for i in range(len(B) - 2, -1, -1):
        X[i] = d_p[i] - c_p[i] * X[i + 1]
    
    
    return X