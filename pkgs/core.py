#!/usr/bin/env python3

from pkgs.riemannian_geometry import diagonal_metric_tensor


def main():
    G = diagonal_metric_tensor(3)

    with open("brew/scalar.txt", "w") as t_file:
        print(G.scalar, file=t_file)


if __name__ == "__main__":
    main()
