from sympy import Function, Matrix, diff, exp, var


class diagonal_metric_tensor:
    def __init__(self, dimension):
        self.n = dimension
        self.x = self.__coords()
        self.g = self.__metric()
        self.h = self.__metric_inv()
        self.Gamma = self.__christoffel_symbol
        self.scalar = self.__scalar()

    def __kronecker_delta(self, i, j):
        return 1 if i == j else 0

    def __coords(self):
        return var(f"X:{self.n}", real=True)

    def __metric(self):
        U = [Function(f"U{i}")(*self.x) for i in range(self.n)]
        G = Matrix(
            [
                [
                    exp(2 * U[j]) * self.__kronecker_delta(i, j)
                    for j in range(self.n)
                ]
                for i in range(self.n)
            ]
        )
        return G

    def __metric_inv(self):
        return self.g**-1

    def __christoffel_symbol(self, i, j, k):
        C = 0
        for l in range(self.n):
            C += (
                0.5
                * (
                    diff(self.g[j, l], self.x[i])
                    + diff(self.g[l, i], self.x[j])
                    - diff(self.g[i, j], self.x[l])
                )
                * self.h[l, k]
            )
        return C.simplify()

    def riemann(self, i, j, k, l):
        R = (
            diff(self.Gamma(i, k, l), self.x[j], 1)
            - diff(self.Gamma(i, j, l), self.x[k], 1)
        )
        for m in range(self.n):
            R += (
                self.Gamma(i, k, m) * self.Gamma(m, j, l)
                - self.Gamma(i, j, m) * self.Gamma(m, k, l)
            )
        return R.simplify()

    def ricci(self, i, j):
        Rc = 0
        for k in range(self.n):
            Rc += self.riemann(i, k, j, k)
        return Rc.simplify()

    def __scalar(self):
        S = 0
        for i in range(self.n):
            for j in range(self.n):
                S += self.ricci(i, j) * self.h[j, i]
        return S.simplify()
