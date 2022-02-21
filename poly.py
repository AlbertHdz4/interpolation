""" Alberto Hernández López"""
""" Interpolacion Polinomial """

import numpy                as np
import sympy                as sym
import matplotlib.pyplot    as plt

import time

from colors             import bcolors
from gauss_elimination  import gauss_elimination


def poly (x_i, f_i) :
    """
    Polinomio de interpolacion
    
    x_i : entradas para la funcion
    f_i : salidas de la funcion 
    returrn : Polinomio de interrpolacion
    """
    n = len(x_i)

    initial_time = time.time()

    # Matriz Vandermonde V
    V = np.zeros((n, n), float)

    for i in range(0, n):
        for j in range(0, n):
            pow     = (n - 1) - j
            V[i,j]  = x_i[i] ** pow


    # Realizamos eliminacion de Gauss - diagonalizacion
    # coefficient = gauss_elimination(V, f_i)
    coefficient = np.linalg.solve(V, f_i)
    # print(coefficient)

    # Usamos Symbol para indicar
    # que x es una variable simbolica
    x = sym.Symbol('x')

    poly = 0

    # Formamos el polinomio
    for i in range(0, n):
        pow         = (n - 1) - i # Potencia para el termino del polinomio
        poly        += coefficient[i] * (x ** pow) # Construimos nuestro polinomio

    """ Con ayuda de sympy se expande el polinomio 
        para su posterior evaluacion, ejemplo de lo que hace
        sympy.expand((x + y)**2) = x**2 +2*x*y + y**2 y
        lambdify es para evaluar todos los puntos de 
        nuestro dominio en el polinomio de Lagrange """
    expanded_interpol_poly = sym.lambdify(x, poly.expand())

    """ Evaluamos nuestra expresion del polinomio 
        en cada punto de nuestras x_i """
    interpol_poly_i = expanded_interpol_poly(x_i)

    final_time  = time.time()
    diff        = final_time - initial_time
    print(f"{bcolors.WARNING}" 'Polinomio de interpolacion: El tiempo que tardo en solucionar el problema con ', n, ' puntos fue de ', diff, ' segundos')

    """ Imprimimos el polinomio """
    # print('Polinomio de interpolación: ')
    # sym.pprint(expanded_interpol_poly)

    """ Grafica - Polinomio de interpolación """
    plt.figure('Polinomio de interpolación')
    plt.title('Polinomio de interpolación')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.plot(x_i, f_i, 'o', color = 'red', label = 'Datos de la muestra')
    plt.plot(x_i, interpol_poly_i, label = poly)
    plt.legend()
    plt.show()

    return poly