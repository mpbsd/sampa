#!/usr/bin/env python3

from sympy import Function, diff, exp, var

DIM = 3

X = var(f"X:{DIM}")
U = Function("U")(*X)


def delta(i, j):
    return 1 if i == j else 0


def g(i, j, *X):
    return exp(2 * U) * delta(i, j)


def h(i, j, *X):
    return exp(-2 * U) * delta(i, j)


def GAMMA(i, j, k, *X):
    G = 0
    for l in range(DIM):
        G += (
            diff(g(j, l, *X), X[i], 1)
            + diff(g(l, i, *X), X[j], 1)
            - diff(g(i, j, *X), X[l], 1)
        ) * h(l, k, *X)
    return G


def RIEM13(i, j, k, l, *X):
    R = (
        diff(GAMMA(i, k, l, *X), X[j], 1)
        - diff(GAMMA(i, j, l, *X), X[k], 1)
    )
    for m in range(DIM):
        R += (
            GAMMA(i, k, m, *X) * GAMMA(m, j, l, *X)
            - GAMMA(i, j, m, *X) * GAMMA(m, k, l, *X)
        )
    return R


def RIEM04(i, j, k, l, *X):
    R = 0
    for m in range(DIM):
        R += RIEM13(i, j, k, m, *X) * g(m, l, *X)
    return R


def RICCI(i, j, *X):
    R = 0
    for k in range(DIM):
        R += RIEM13(i, k, j, k, *X)
    return R


def SCAL(*X):
    S = 0
    for i in range(DIM):
        for j in range(DIM):
            S += RICCI(i, j, *X) * h(j, i, *X)
    return S


def main():
    print(SCAL(*X))


if __name__ == "__main__":
    main()
