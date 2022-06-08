import math

import sympy

"""
 * Authors: Orel Dandeker 
"""

from math import e
import numpy as np
import sympy as sp
import math
from sympy.utilities.lambdify import lambdify

x = sp.symbols('x')


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[90m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    # Background colors:
    GREYBG = '\033[100m'
    REDBG = '\033[101m'
    GREENBG = '\033[102m'
    YELLOWBG = '\033[103m'
    BLUEBG = '\033[104m'
    PINKBG = '\033[105m'
    CYANBG = '\033[106m'


def TrapezoidalRule(f, n, a, b, tf):
    """
    rapezoidal Rule is a rule that evaluates the area under the curves by dividing the total area
    into smaller trapezoids rather than using rectangles
    :param f: The desired integral function
    :param n: The division number
    :param a: Lower bound
    :param b: Upper bound
    :param tf: Variable to decide whether to perform Error evaluation
    :return: The result of the integral calculation
    """
    h = (b - a) / n
    if tf:
        print(bcolors.FAIL, "Error evaluation En = ", round(TrapezError(func(), b, a, h), 6), bcolors.ENDC)
    integral = 0.5 * (f(a) * f(b))
    for i in range(n):
        integral += f(a + h * i)
    integral *= h
    return integral


def SimpsonRule(f, n, a, b):
    """
    Simpson’s Rule is a numerical method that approximates the value of a definite integral by using quadratic
     functions Simpson’s Rule is based on the fact that given three points,
    we can find the equation of a quadratic through those points (by Lagrange's interpolation)
    :param f: The desired integral function
    :param n: The division number(must be even)
    :param a: Lower bound
    :param b: Upper bound
    :return: The result of the integral calculation
    """
    if n % 2 != 0:
        return 0, False
    h = (b - a) / n
    print(bcolors.FAIL, "Error evaluation En = ", round(SimpsonError(func(), b, a, h), 6), bcolors.ENDC)
    integral = f(a) + f(b)
    for i in range(n):
        k = a + i * h
        if i % 2 == 0:
            integral += 2 * f(k)
        else:
            integral += 4 * f(k)
    integral *= (h / 3)
    return integral, True


def RombergsMethod(f, n, a, b):
    """
    Romberg integration is an extrapolation technique which allows us to take a sequence
    approximate solutions to an integral and calculate a better approximation.
    This technique assumes that the function we are integrating is sufficiently differentiable
    :param f: The desired integral function
    :param n: The division number
    :param a: Lower bound
    :param b: Upper bound
    :return: The result of the integral calculation
    """
    matrix = [[0 for i in range(n)] for j in range(n)]
    for k in range(0, n):
        # Using the trapezoidal method
        matrix[k][0] = TrapezoidalRule(f, 2 ** k, a, b, False)
        # Romberg recursive formula Using values that have already been calculated
        for j in range(0, k):
            matrix[k][j + 1] = (4 ** (j + 1) * matrix[k][j] - matrix[k - 1][j]) / (4 ** (j + 1) - 1)
            print("R[{0}][{1}] = ".format(k, j + 1), round(matrix[k][j + 1], 6))
    return matrix


def func():
    return (x * math.e ** (-x) + sp.ln(x ** 2)) * (2 * x ** 3 + 2 * x ** 2 - 3 * x - 5)


def f(val):
    return lambdify(x, func())(val)


def TrapezError(func, b, a, h):
    """
    The trapezoidal rule is a method for approximating definite integrals of functions.
    The error in approximating the integral of a twice-differentiable function by the trapezoidal rule
    is proportional to the second derivative of the function at some point in the interval.
    :param func: The desired integral function
    :param b: Upper bound
    :param a: Lower bound
    :param h: The division
    :return: The error value
    """
    xsi = (-1) * (math.pi / 2)
    print("ƒ(x): ", func)
    f2 = sp.diff(func, x, 2)
    print("ƒ'': ", f2)
    diff_2 = lambdify(x, f2)
    print("ƒ''(", xsi, ") =", diff_2(xsi))
    return h ** 2 / 12 * (b - a) * diff_2(xsi)


def SimpsonError(func, b, a, h):
    """
    The Simpson rule is a method for approximating definite integrals of functions.
    The error in approximating the integral of a four-differentiable function by the trapezoidal rule
    is proportional to the second derivative of the function at some point in the interval.
    :param func: The desired integral function
    :param b: Upper bound
    :param a: Lower bound
    :param h: The division
    :return: The error value
    """
    xsi = 1
    print("ƒ(x): ", func)
    f2 = sp.diff(func, x, 4)
    print("ƒ⁴: ", f2)
    diff_4 = lambdify(x, f2)
    print("ƒ⁴(", xsi, ") =", diff_4(xsi))

    return (math.pow(h, 4) / 180) * (b - a) * diff_4(xsi)


def MainFunction():
    n = 2
    print(bcolors.BOLD, "Division into sections n =", n, bcolors.ENDC)
    print(bcolors.OKBLUE, "Numerical Integration of definite integral in range [0,𝛑] ∫= SIN(X)", bcolors.ENDC)
    choice = int(input(
        "Which method do you want? \n\t1.The Trapezoidal Rule \n\t2.Simpson’s Rule\n\t3.Romberg's method\n"))
    if choice == 1:
        print(bcolors.OKBLUE, "I = ", round(TrapezoidalRule(f, n, 0, math.pi, True), 6), bcolors.ENDC)
    elif choice == 2:
        res = SimpsonRule(f, n, 0.5, 1)
        if res[1]:
            print(bcolors.OKBLUE, "I = ", round(res[0], 6), bcolors.ENDC)
        else:
            print(bcolors.FAIL, "n must be even !", bcolors.ENDC)
    elif choice == 3:
        print(bcolors.OKBLUE, "I = ", round(RombergsMethod(f, n, 0, math.pi)[n - 1][n - 1], 6), bcolors.ENDC)
    else:
        print(bcolors.FAIL, "Invalid input", bcolors.ENDC)


MainFunction()
