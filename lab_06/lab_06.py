import math

import numpy 

# CONST

Eps = 1e-4

def f(t): return math.exp(-t * t / 2)

def GetPolinomLegendre(n, x):
    p = []
    p.append(1)
    p.append(x)

    i = 2

    while (i <= n):
        tmp = (2 * i - 1) * x * p[i - 1] - (i -1) * p[i - 2]
        tmp /= i
        p.append(tmp)

        i += 1

    return p

def GetDiveration(n, p, x): return n / (1 - x * x) * (p[n -1] - x * p[n])

def RootsLegendre(n):
    px =[]
    dpx = []

    x = [math.cos(math.pi * (4*i-1) / (4*n+2)) for i in range(n, 0, -1)]

    for i in range(0, n):
        p = []
        dp = []

        while (True):
            p = GetPolinomLegendre(n, x[i])
            dp = GetDiveration(n, p, x[i])
            dx = p[n] / dp

            x[i] -= dx

            if math.fabs(dx) < Eps:
                break

        px.append(p)
        dpx.append(dp)

    return x, px, dpx

def GaussFormulaCoefs(x, dp, n):

    def GetAlpha(i):
        return 2 / (1 - x[i]*x[i]) / (dp[i]**2)

    Alpha = [GetAlpha(i) for i in range(0, n)]

    return Alpha

def GaussCoefsB(x, n):
    z = []

    for i in range(0, n):
        if (i % 2 == 0):
            z.append(2 / (i + 1))
        else:
            z.append(0)

    matrix = []

    matrix.append([1 for i in range(0, n)])

    for i in range(1, n):
        matrix.append([])
        for j in range(0, n):
            matrix.append(matrix[i -1][j] * x[j])

    result = numpy.linalg.solve(matrix, z)

    return result

def F(x, Alpha, t, weights):
    result = 0

    n = len(t)

    for i in range(0, n):
        result += weights[i] * f((x / 2) * (t[i] + 1))

    result *= (x / 2)

    return result - Alpha

def FindIntegrationLimit(a, b, Alpha, t, weights):
    if (F(a, Alpha, t, weights) > 0):
        a, b = b, a

    tmp = 0

    while (True):
        tmp = (a + b) / 2
        Ftmp = F(tmp, Alpha, t, weights)

        if (math.abs((b - tmp) / b) < Eps):
            break

        if Ftmp < 0: a = tmp
        else: b = tmp

    return tmp

def Main():
    n = int(input('n = '))
    Alpha = float(input('a = '))

    Alpha = Alpha * math.sqrt(2 * math.pi)

    x, px, dpx = RootsLegendre(n)

    #print("Roots of Legendre's polynome >> ", x, sep='\n')

    #print("Weights coefficients >> ", GaussFormulaCoefs(x, dpx, n), sep='\n')

    b = GaussCoefsB(x, n)

    result = FindIntegrationLimit(0, 5, Alpha, x, b)
    print('x = ', result)
    print(F(result, Alpha, x, b))

if __name__ == "__main__":
    Main()
