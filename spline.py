"""" Alberto Hernández López """
"""" Interpolacion de splines lineal 
    y cubico para funciones 1D """

import matplotlib.pyplot    as plt
import numpy                as np
import sympy                as sym
import bisect
import math
import time

from gauss_elimination  import gauss_elimination
from colors             import bcolors


# Spline lineal para la interpolacion de funciones 1D
def linear_spline (x_i, f_i) :
    """
    Interpola funciones 1D usando un spline lineal
    x_i : Coordenadas x a ser evaluadas en f(x) y a ser interpolados 
    f_i : Valores o salidas de la funcion

    Se deben de sumar polinomio de la forma, 
    Pi(x) = ax + b
    """
    initial_time = time.time()

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
    
    final_time  = time.time()
    diff        = final_time - initial_time
    print(f"{bcolors.WARNING}" 'Spline lineal: El tiempo que tardo en solucionar el problema con ', n, ' puntos fue de ', diff, ' segundos')

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
    

def get_cubic_spline (x_i, f_i) :
    """
    Interpola funciones 1D usando un spline cubico
    x_i : Coordenadas x a ser evaluadas en f(x) y a ser interpolados 
    f_i : Valores o salidas de la funcion

    Se deben de sumar polinomio de la forma, 
    Pi(x) = ax**3 + bx**2 + cx + d
    """
    n = len(x_i)

    # Creamos la matriz h
    h = np.zeros(n - 1, float)

    # Formamos la matriz h
    for j in range(0, n - 1) :
        h[j] = x_i[j+1] - x_i[j]
    
    # Declaramos las matrices que resolveremos por gauss
    A = np.zeros((n - 2, n - 2), float)
    B = np.zeros(n - 2, float)
    S = np.zeros(n, float)
    
    # Asignamos los primeros elementos
    A[0, 0] = 2 * (h[0] + h[1])
    A[0, 1] = h[1]
    B[0]    = 6 * ((f_i[2] - f_i[1]) / h[1] - (f_i[1] - f_i[0]) / h[0])

    # Empezamos a formar las matrices
    for i in range(1, n-3):
        A[i, i - 1] = h[i]
        A[i, i]     = 2 * (h[i] + h[i + 1])
        A[i, i + 1] = h[i + 1]
        factor_2_1  = (f_i[i + 2] - f_i[i + 1]) / h[i + 1]
        factor_1_0  = (f_i[i + 1] - f_i[i]) / h[i]
        B[i]        = 6 * (factor_2_1 - factor_1_0)
        
    A[n - 3, n - 4] = h[n - 3]
    A[n - 3, n - 3] = 2 * (h[n - 3] + h[n - 2])
    factor_1_2  = (f_i[n - 1] - f_i[n - 2]) / h[n - 2]
    factor_2_3  = (f_i[n - 2] - f_i[n - 3]) / h[n - 3]
    B[n - 3]    = 6 * (factor_1_2 - factor_2_3)
    
    solutions = gauss_elimination(A, B)
    
    # Soluciones con linalg solo para comprobar que los resultados de
    # gauss_elimination sean correctos
    solutions = np.linalg.solve(A, B) 

    for j in range(1, n - 1) :
        S[j] = solutions[j - 1]

    S[0]        = 0
    S[n - 1]    = 0
    

    a_c = np.zeros(n-1, dtype = float)
    b = np.zeros(n-1, dtype = float)
    c = np.zeros(n-1, dtype = float)
    d = np.zeros(n-1, dtype = float)

    for j in range(0, n - 1) :
        a_c[j]  = (S[j + 1] - S[j]) / (6 * h[j])
        b[j]    = S[j] / 2
        factor_1_0 = (f_i[j + 1] - f_i[j]) / h[j]
        c[j] = factor_1_0 - (2 * h[j] * S[j] + h[j] * S[j + 1]) / 6
        d[j] = f_i[j]
    
    """ Indicamos que el x sera simbolica y la evaluaremos posteriormente"""
    x = sym.Symbol('x')
    
    all_cubic_polys = []

    for j in range(0,n-1,1):
        cubic_poly = a_c[j] * (x - x_i[j]) ** 3 + b[j] * (x - x_i[j]) ** 2
        cubic_poly = cubic_poly + c[j] * (x - x_i[j]) + d[j]
        
        cubic_poly = cubic_poly.expand()
        all_cubic_polys.append(cubic_poly)
            
    """Graficamos los polinomios cubicos con sus respectivos dominios """
    x_poly_i    = np.array([])
    y_poly_i    = np.array([])
    step        = 1

    while not(step >= n):
        x_temp_i = x_i[step - 1]
        y_temp_i = x_i[step]
        x_temp   = np.linspace(x_temp_i, y_temp_i, n)

        # evalua polinomio del tramo
        cubic_poly_i = all_cubic_polys[step - 1]
        cubic_poly_t = sym.lambdify('x', cubic_poly_i) 
        f_i_t = cubic_poly_t(x_temp)

        # vectores de trazador en x,y
        x_poly_i = np.concatenate((x_poly_i, x_temp))
        y_poly_i = np.concatenate((y_poly_i, f_i_t))
        step = step + 1

    """ Grafiamos los polinomios cubicos """
    plt.plot(x_i, f_i, 'o', color = 'red', label = 'Puntos a ajustarr')
    plt.plot(x_poly_i, y_poly_i, label = 'Spline cubico')
    plt.title('Spline cubico')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.show()

    return(all_cubic_polys)