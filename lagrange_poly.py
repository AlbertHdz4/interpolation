""" Alberto Hernandez Lopez """
"""" Polinomio de Lagrange """

import matplotlib.pyplot    as plt
import numpy                as np
import sympy                as sym
import time

from colors             import bcolors


def lagrange_poly (x_i, f_i) :
    """
    Polinomio de Lagrange
    
    x_i     : entradas de la funcion, puntos x
    f_i     : salidas de la funcion a interpolar
    return  : regresamos el polinomio de Lagrange
    """
    initial_time = time.time()
    
    n = len(x_i)
    x = sym.Symbol('x') # Indicamos que x sera nuestra variable a simbolica a evaluar

    # Polinomio de Lagrange
    lagrange_poly = 0

    # Empezamos a construir los terminos del polinomio de Lagrange
    for i in range(0, n):
        numerator   = 1
        denominator = 1
        for j in range(0,n):
            if(i!=j):
                numerator = numerator * (x - x_i[j])
                denominator = denominator * (x_i[i] - x_i[j])
            lagrangian_term = (numerator/denominator) * f_i[i]
        lagrange_poly = lagrange_poly + lagrangian_term


    """ Con ayuda de sympy se expande el polinomio 
    para su posterior evaluacion, ejemplo de lo que hace
    sympy.expand((x + y)**2) = x**2 +2*x*y + y**2 y
    lambdify es para evaluar todos los puntos de 
    nuestro dominio en el polinomio de Lagrange """
    simbolic_poly = sym.lambdify(x, lagrange_poly)


    """ Evaluamos los x_i en el polinomio"""
    x_min_new   = np.min(x_i)
    x_max_new   = np.max(x_i)
    p_x_i       = np.linspace(x_min_new, x_max_new, n)
    lagrange_poly_f_i = simbolic_poly(p_x_i)


    """ Obtenemos los tiempos finales que tardo el algoritmo """
    final_time  = time.time()
    diff        = final_time - initial_time
    print(f"{bcolors.WARNING}"'Polinomio de Lagrange: El tiempo que tardo en solucionar el problema con ', n, ' puntos fue de ', diff, ' segundos')


    """ Graficamos el polinomio de Lagrange """
    plt.figure('Interpolación Lagrange')
    plt.title('Interpolación Lagrange')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.plot(x_i, f_i, 'o', color = 'red', label = 'Puntos de la muestra')
    plt.plot(p_x_i, lagrange_poly_f_i, label = 'Polinomio de Lagrange')
    plt.legend()
    plt.show()

    return lagrange_poly