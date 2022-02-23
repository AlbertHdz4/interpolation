""" Alberto Hernández López"""
""" Interpolacion """

# Imports
import numpy                as np # Para el manejo de matrices
import matplotlib.pyplot    as plt # Para graficar
import random               as rd # Para generar numeros aleatorios
import time


# Sympy es utilizado para expresar de manera 
# más rápida y sencilla los polinomios y sus 
# evaluaciones https://docs.sympy.org/latest/modules/utilities/lambdify.html
import sympy                as sym 

from lagrange_poly          import * # Polinomio de Lagrange
from poly                   import * # Polinomio de interpolacion
from spline                 import * # Spline lineal y cubico
from colors                 import * # Para dar estilo en terminal


## Para generar un rango aleatorio
def space_gen(x_min, x_max, samples):
    x = np.linspace(x_min, x_max, samples, float)
    for i in range(len(x) - 1) :
        rd_number   = rd.uniform(0, 1)
        x[i] *= rd_number 

    x.sort()
    x[0]            = x_min
    x[len(x) - 1]   = x_max
    return x


if __name__ == '__main__':
    print(f"{bcolors.OKCYAN}******************** START ********************{bcolors.ENDC}")

    samples = 16 # 16, 32, 64 # Definimos los numeros de puntos deseados


    print(f"{bcolors.OKGREEN}******************** PRIMERA PARTE A) ********************{bcolors.ENDC}")
    """ Interpolación con espaciamiento 
        regular e irregular de 16, 32, 64 puntos """
    # (x_i, step) = np.linspace(0, 360, samples, retstep = True) # Espaciamiento regular de 0 - 2*Pi
    # (x_i, step) = np.linspace(0, 1, samples, dtype = float, retstep = True) # Espaciamiento regular de 0 - 1
    x_i = space_gen(0, 2 * np.pi, samples) # Espaciamiento aleatorio de 0 - 2*Pi 
    # x_i = space_gen(0, 1, samples) # Espaciamiento aleatorio de 0 - 1


    """ Funcion f(x) """
    f_i = 2 * np.cos(x_i) + np.sin(x_i) + np.sqrt(x_i)
    # f_i = 2 * np.cos(2 * np.pi * x_i) + np.sin(2 * np.pi * x_i) + np.sqrt(2 * np.pi * x_i)

    
    """ Grafica - puntos a ajustar ------------------------------------------"""
    """" Solo vemos los puntos"""
    plt.figure("Puntos a interpolar")
    plt.title("Puntos a ajustar")
    plt.xlabel('x')
    plt.xlabel('y')
    plt.plot(x_i, f_i, 'x', color = "purple") 
    plt.show()
    """ ---------------------------------------------------------------------- """

    """ Si deseas probar, tienes que descomentar un dominio, un rango y el metodo
        de interpolacion o spline que desees, todo esta seccionado por metodo """

    """ Polinomio de interpolacion -------------- """
    dom_interpol_poly       = np.linspace(0, 360, samples) # Dominio regular de 0 - 2 * Pi
    # dom_interpol_poly     = np.linspace(0, 1, samples, dtype = float) # Espaciamiento regular de 0 - 1
    # dom_interpol_poly     = space_gen(0, 360, samples) # Dominio con espaciamiento aleatorio de 0 - 2 * Pi
    # dom_interpol_poly     = space_gen(0, 1, samples) # Dominio con espaciamiento aleatorio de 0 - 1
    range_interpol_poly     = 2 * np.cos(dom_interpol_poly) + np.sin(dom_interpol_poly) + np.sqrt(dom_interpol_poly) # Funcion a aproximar
    # range_interpol_poly   = 2 * np.cos(2 * np.pi * dom_interpol_poly) + np.sin(2 * np.pi * dom_interpol_poly) + np.sqrt(2 * np.pi * dom_interpol_poly) # Otra funcion a aproximar
    inter_poly          = poly(dom_interpol_poly, range_interpol_poly)   
    """ ---------------------------------------------------------------------- """


    """ Polinomio de Lagrange ------------------------------------------------ """
    # dom_lagrange  = np.linspace(0, 360, samples) # Dominio regular de 0 - 2*Pi
    # dom_lagrange  = np.linspace(0, 1, samples, dtype = float) # Espaciamiento regular de 0 - 1
    # dom_lagrange  = space_gen(0, 360, samples) # Dominio con espaciamiento aleatorio de 0 - 2*Pi
    # dom_lagrange  = space_gen(0, 1, samples) # Dominio con espaciamiento aleatorio de 0 - 1
    # range_lagrange  = 2 * np.cos(dom_lagrange) + np.sin(dom_lagrange) + np.sqrt(dom_lagrange) # Funcion a aproximar
    # range_lagrange  = 2 * np.cos(2 * np.pi * dom_lagrange) + np.sin(2 * np.pi * dom_lagrange) + np.sqrt(2 * np.pi * dom_lagrange) # Otra funcion a aproximar
    # lagrangian_polynomial = lagrange_poly(dom_lagrange, range_lagrange)
    """ ---------------------------------------------------------------------- """


    print(f"{bcolors.OKGREEN}\n\n******************** SEGUNDA PARTE ********************{bcolors.ENDC}")
    """ Spline lineal ------------------------------------------------ """
    # dom_linear_spline     = np.linspace(0, 360, samples) # Dominio regular de 0 - 2*Pi
    # dom_linear_spline     = np.linspace(0, 1, samples, dtype = float) # Espaciamiento regular de 0 - 1
    dom_linear_spline     = space_gen(0, 360, samples) # Dominio con espaciamiento aleatorio de 0 - 2*Pi
    # dom_linear_spline     = space_gen(0, 1, samples) # Dominio con espaciamiento aleatorio de 0 - 1
    range_linear_spline   = 2 * np.cos(dom_linear_spline) + np.sin(dom_linear_spline) + np.sqrt(dom_linear_spline) # Funcion a aproximar
    # range_linear_spline   = 2 * np.cos(2 * np.pi * dom_linear_spline) + np.sin(2 * np.pi * dom_linear_spline) + np.sqrt(2 * np.pi * dom_linear_spline) # Otra funcion a aproximar
    linear_poly = linear_spline(dom_linear_spline, range_linear_spline)
    """ ---------------------------------------------------------------------- """


    """ Spline cubico ------------------------------------------------ """
    # dom_cubic_spline  = np.linspace(0, 360, samples) # Dominio regular de 0 - 2*Pi
    # dom_cubic_spline  = np.linspace(0, 1, samples, dtype = float) # Espaciamiento regular de 0 - 1
    # dom_cubic_spline  = space_gen(0, 360, samples) # Dominio con espaciamiento aleatorio de 0 - 2*Pi
    # dom_cubic_spline    = space_gen(0, 1, samples) # Dominio con espaciamiento aleatorio de 0 - 1
    # range_cubic_spline  = 2 * np.cos(dom_cubic_spline) + np.sin(dom_cubic_spline) + np.sqrt(dom_cubic_spline) # Funcion a aproximar
    # range_cubic_spline  = 2 * np.cos(2 * np.pi * dom_cubic_spline) + np.sin(2 * np.pi * dom_cubic_spline) + np.sqrt(2 * np.pi * dom_cubic_spline) # Otra funcion a aproximar
    # cubic_spline        = get_cubic_spline(dom_cubic_spline, range_cubic_spline)
    """ ---------------------------------------------------------------------- """
    ## ----------------------------------------------------------------------
    print(f"{bcolors.OKCYAN}\n******************** END ********************{bcolors.ENDC}")
