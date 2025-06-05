#!/usr/bin/env python3

from sympy import Function, diff, exp, latex, var

DIM = 4

X = var(f"X:{DIM}")
U = Function("U")(*X)


def delta(i, j):
    return 1 if i == j else 0


def g(i, j, *X):
    """
    g_{ij} = g(\partial_{i},\partial_{j})
    """
    return exp(2 * U) * delta(i, j)


def h(i, j, *X):
    """
    \sum_{k=1}^{DIM}g_{ik}h_{kj} = \delta_{ij}
    """
    return exp(-2 * U) * delta(i, j)


def gamma(i, j, k, *X):
    """
    \\nabla_{\partial_{i}}\partial_{j}
    =
    \sum_{k=1}^{DIM}\Gamma_{ij}^{k}\partial_{k}
    """
    G = 0
    for l in range(DIM):
        G += (
            diff(g(j, l, *X), X[i], 1)
            + diff(g(l, i, *X), X[j], 1)
            - diff(g(i, j, *X), X[l], 1)
        ) * h(l, k, *X)
    return G


def riem13(i, j, k, l, *X):
    """
    Returns the R_{ijk}^{l} in:
    R(\partial_{i},\partial{j})\partial_{k}
    =
    \sum_{l=1}^{DIM}R_{ijk}^{l}\partial_{l}
    """
    R = (
        diff(gamma(i, k, l, *X), X[j], 1)
        - diff(gamma(i, j, l, *X), X[k], 1)
    )
    for m in range(DIM):
        R += (
            gamma(i, k, m, *X) * gamma(m, j, l, *X)
            - gamma(i, j, m, *X) * gamma(m, k, l, *X)
        )
    return R


def riem04(i, j, k, l, *X):
    """
    Returns g(R(\partial_{i},\partial_{j},\partial_{k}),\partial_{l}):
    """
    R = 0
    for m in range(DIM):
        R += riem13(i, j, k, m, *X) * g(m, l, *X)
    return R


def ricci(i, j, *X):
    R = 0
    for k in range(DIM):
        R += riem13(i, k, j, k, *X)
    return R


def scal(*X):
    S = 0
    for i in range(DIM):
        for j in range(DIM):
            S += ricci(i, j, *X) * h(j, i, *X)
    return S


def main():
    print(scal(*X))


if __name__ == "__main__":
    main()
