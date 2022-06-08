"""
 * Authors: Orel Dandeker
"""

import math

import sympy
import sympy as sp
from sympy.utilities.lambdify import lambdify


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[90m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def Run(choice, pol, pol_prime):
    if choice == 1:
        print(bcolors.OKBLUE, "Bisection method on: ", pol, bcolors.ENDC)
        Call_Method(F, 1)
        ranges_with_roots.clear()
        print(bcolors.OKBLUE, "\nFinding roots on: ", pol_prime, bcolors.ENDC)
        find_roots_ranges(Range[0], Range[1], F_prime)
        PrintingRoots()
        print(bcolors.OKBLUE, "\nBisection method on: ", pol_prime, bcolors.ENDC)
        Call_Method(F_prime, 1)
        if F(0) == 0:
            print(bcolors.UNDERLINE, "0 is also root of: ", pol, "( F(0)=0 ) ", bcolors.ENDC)
        ranges_with_roots.clear()
    elif choice == 2:
        print(bcolors.OKBLUE, "Newton-Raphson method on: ", pol, bcolors.ENDC)
        Call_Method(F, 2)
        if F(0) == 0:
            print(bcolors.UNDERLINE, "0 is also root of: ", pol, "( F(0)=0 ) ", bcolors.ENDC)
    else:
        print(bcolors.OKBLUE, "Secant method on: ", pol, bcolors.ENDC)
        Call_Method(F, 3)
        if F(0) == 0:
            print(bcolors.UNDERLINE, "0 is also root of: ", pol, "( F(0)=0 ) ", bcolors.ENDC)


def PrintingRoots():
    if len(ranges_with_roots) != 0:
        print(bcolors.OKBLUE, "The ranges with the roots found are:", bcolors.ENDC)
        for i in range(int(len(ranges_with_roots) / 2)):
            if i == 0:
                print((ranges_with_roots[0], ranges_with_roots[1]))
            else:
                print((ranges_with_roots[i * 2], ranges_with_roots[(i * 2) + 1]))
    elif len(ranges_with_roots) == 0 and len(roots) != 0:
        print(bcolors.HEADER, "Sign changing not found in the specific selected range !\nThe only root/s found: ",
              roots, bcolors.HEADER)
        roots.clear()
    else:
        print(bcolors.HEADER, "No root found in the specific selected range ! ", bcolors.ENDC)


def calculate_Bisection_Method_Eror(eror, a, b):
    x = math.log((eror) / (b - a))
    x = math.ceil((-x) / math.log(2))
    return x


def Bisection_Method(polinom, start, end, MaxIterations, eps=0.00001):
    A = round(start, 2)
    B = end
    counter = 0
    c = 0
    while (B - A) > eps:
        counter += 1
        if counter == MaxIterations:
            if round((A + B) / 2, 5) == c and round(abs(F(c)), 2) == 0:
                print(bcolors.HEADER, "The root found is: ", c, bcolors.ENDC)
                print("\n")
                return
            elif round((A + B) / 2, 5) == c and round(abs(F(c)), 2) != 0:
                print(bcolors.HEADER, "The root found is:", c, "but this root is not suitable ! F(", c,
                      ") is not equal to 0", bcolors.ENDC)
                print("\n")
                return
            else:
                print(bcolors.HEADER,
                      "The BisectionMethod is not suitable ! We reached to max iterations and no root found",
                      bcolors.ENDC)
                return
        c = (A + B) / 2
        c = round(c, 5)
        print("iteration number ", counter, "\t", "A=", A, "\t", "B=", B, "\t", "C=", c)
        if polinom(A) * polinom(c) < 0:
            B = c
        else:
            A = c

    if abs(round(F(c), 2)) == 0:
        print(bcolors.HEADER, "The root found is: ", c, bcolors.ENDC)
        print("\n")
    else:
        print(bcolors.FAIL, "The root found is:", c, "but this root is not suitable ! F(", c, ") is not equal to 0",
              bcolors.ENDC)
        print("\n")


def find_roots_ranges(a, b, F):
    i = 1
    while a <= b:
        a = round(a, 2)
        if round(F(a), 3) == 0:
            roots.append(a)
        print("iteration number ", i, "\t", "x=", a, "\t", "F(x)=", F(a))
        if F(a) * F(a - 0.1) < 0:
            print(bcolors.UNDERLINE, "sign changing has been found between iterations ", i, "-", i - 1, bcolors.ENDC)
            ranges_with_roots.append(round(a - 0.1, 2))
            ranges_with_roots.append(round(a, 2))
        i += 1
        a += 0.1
    print("\n")


def Call_Method(F, methodNum):
    for i in range(int(len(ranges_with_roots) / 2)):
        if methodNum == 1:
            if i == 0:
                Bisection_Method(F, ranges_with_roots[0], ranges_with_roots[1], maxIterations)
            else:
                Bisection_Method(F, ranges_with_roots[i * 2], ranges_with_roots[(i * 2) + 1], maxIterations)
        elif methodNum == 2:
            if i == 0:
                NewtonRaphson_Method(F, ranges_with_roots[0], ranges_with_roots[1], maxIterations)
            else:
                NewtonRaphson_Method(F, ranges_with_roots[i * 2], ranges_with_roots[(i * 2) + 1], maxIterations)
        else:
            if i == 0:
                Secant_Method(F, ranges_with_roots[0], ranges_with_roots[1], maxIterations)
            else:
                Secant_Method(F, ranges_with_roots[i * 2], ranges_with_roots[(i * 2) + 1], maxIterations)


def Secant_Method(polinom, firstGuess, secondGuess, maxIterations, eps=0.00001):
    Xi = firstGuess
    Xi_1 = secondGuess
    counter = 0
    print("Start with: \t Xi(firstGuess)=", Xi, "\t", "Xi_1(secondGuess)=", Xi_1, "\t", "F(Xi)=", round(polinom(Xi), 6))
    while abs(Xi_1 - Xi) > eps:
        counter += 1
        if counter == maxIterations:
            if round((Xi * polinom(Xi_1) - Xi_1 * polinom(Xi)) / (
                    polinom(Xi_1) - polinom(Xi)), 5) == Xi_1 and round(abs(F(Xi_1)), 2) == 0:
                print(bcolors.HEADER, "The root found is: ", Xi_1, bcolors.ENDC)
                print("\n")
                return
            else:
                print(bcolors.FAIL, "The SecantMethod is not suitable ! We reached to max iterations and no root found",
                      bcolors.ENDC)
                return
        tmp = Xi_1
        Xi_1 = round((Xi * polinom(Xi_1) - Xi_1 * polinom(Xi)) / (
                polinom(Xi_1) - polinom(Xi)), 5)
        Xi = tmp
        print("iteration number ", counter, "\t", "Xi=", Xi, "\t", "Xi_1=", Xi_1, "\t", "F(Xi)=", round(polinom(Xi), 6))
    print(bcolors.HEADER, "The root found is: ", round(Xi_1, 5), bcolors.ENDC)
    print("\n")


def NewtonRaphson_Method(polinom, start, end, MaxIterations, eps=0.00001):
    counter = 0
    Xr = round((start + end) / 2, 2)
    Xr_1 = Xr - (polinom(Xr) / F_prime(Xr))
    print("Start with :\t Xr=", round(Xr, 6), "\t", "F(x)=", round(polinom(Xr), 6), "\t", "F'(Xr)=",
          round(F_prime(Xr), 6))
    while abs(Xr_1 - Xr) > eps:
        counter += 1
        if counter == maxIterations:
            if round(Xr_1 - (polinom(Xr_1) / F_prime(Xr_1)), 5) == Xr_1 and round(abs(F(Xr_1)), 2) == 0:
                print(bcolors.HEADER, "The root found is: ", Xr_1, bcolors.ENDC)
                print("\n")
                return
            else:
                print(bcolors.FAIL,
                      "The Newton-Raphson Method is not suitable ! We reached to max iterations and no root found",
                      bcolors.ENDC)
                return
        print("iteration number ", counter, "\t", "Xr=", round(Xr_1, 7), "\t", "F(x)=", round(polinom(Xr_1), 6), "\t",
              "F'(Xr)=", round(F_prime(Xr_1), 6))
        tmp = Xr_1
        Xr_1 = Xr_1 - (polinom(Xr_1) / F_prime(Xr_1))
        Xr = tmp
    print(bcolors.HEADER, "The root is: ", round(Xr, 5), bcolors.ENDC)


# ------------------------------Main---------------------------------
roots = []
Range = [0, 1.5]
x = sp.symbols('x')
F = (x * pow(2.718, -x) + sympy.log(x ** 2)) * (2 * x ** 3 + 2 * x ** 2 - 3 * x - 5)
F_prime = sp.diff(F, x)
print(bcolors.OKBLUE, "\nFinding roots on : ", F, bcolors.ENDC)
Fcopy = F
Fprimecopy = F_prime
F = lambdify(x, F)
F_prime = lambdify(x, F_prime)
ranges_with_roots = []
maxIterations = calculate_Bisection_Method_Eror(0.00001, Range[0], Range[1])
find_roots_ranges(Range[0], Range[1], F)
PrintingRoots()
if len(ranges_with_roots) != 0:
    choice = int(input("Choose your method : 1-BisectionMethod\t2-NewtonRaphsonMethod\t3-SecantMethod\n"))
    Run(choice, Fcopy, Fprimecopy)
