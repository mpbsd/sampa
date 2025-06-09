#!/usr/bin/env python3

from sympy import Function, diff, exp, var

DIM = 3

X = var(f"X:{DIM}")
U = [Function(f"U{i}")(*X) for i in range(DIM)]


def delta(i, j):
    return 1 if i == j else 0


def g(i, j):
    return exp(2 * U[j]) * delta(i, j)


def h(i, j):
    return exp(-2 * U[j]) * delta(i, j)


def gamma(i, j, k):
    G = 0
    for l in range(DIM):
        G += (
            diff(g(j, l), X[i], 1)
            + diff(g(l, i), X[j], 1)
            - diff(g(i, j), X[l], 1)
        ) * h(l, k)
    return G


def riem13(i, j, k, l):
    R = (
        diff(gamma(i, k, l), X[j], 1)
        - diff(gamma(i, j, l), X[k], 1)
    )
    for m in range(DIM):
        R += (
            gamma(i, k, m) * gamma(m, j, l)
            - gamma(i, j, m) * gamma(m, k, l)
        )
    return R


def riem04(i, j, k, l):
    R = 0
    for m in range(DIM):
        R += riem13(i, j, k, m) * g(m, l)
    return R


def ricci(i, j):
    R = 0
    for k in range(DIM):
        R += riem13(i, k, j, k)
    return R


def scal():
    S = 0
    for i in range(DIM):
        for j in range(DIM):
            S += ricci(i, j) * h(j, i)
    return S


def main():
    with open("brew/scal.txt", "w") as t_file:
        print(scal(), file=t_file)


if __name__ == "__main__":
    main()
