from __future__ import annotations
from collections import defaultdict
from math import isqrt
from typing import overload


class Poly:
    """
    Poly is a simple polynomial that stores a sum of a*x^b*y^c terms.
    """

    def __init__(
            self,
            constant: int = 0, x_pow: int = 0, y_pow: int = 0) -> None:
        self.coeff: defaultdict[int, defaultdict[int, int]] = defaultdict(
            lambda: defaultdict(int))
        self.set_coeff(constant, x_pow, y_pow)

    def all_coeffs(self) -> list[tuple[int, int, int]]:
        return [
            (coeff, x, y)
            for x in self.coeff.keys()
            for y, coeff in self.coeff[x].items()
            if coeff != 0
        ]

    def get_coeff(self, x_pow: int, y_pow: int) -> int:
        if len(self.coeff) <= x_pow or len(self.coeff[x_pow]) <= y_pow:
            return 0
        return self.coeff[x_pow][y_pow]

    def set_coeff(self, coeff: int, x_pow: int, y_pow: int) -> None:
        if coeff == 0:
            if x_pow not in self.coeff or y_pow not in self.coeff[x_pow]:
                return
            # Delete to save a bit of memory.
            del self.coeff[x_pow][y_pow]
            if len(self.coeff[x_pow]) == 0:
                del self.coeff[x_pow]
            return
        self.coeff[x_pow][y_pow] = coeff

    def highest_x_power(self) -> int:
        return max(self.coeff.keys())

    def highest_y_power(self) -> int:
        return max(y_pow for x_pow in self.coeff.keys()
                   for y_pow in self.coeff[x_pow].keys())

    def __add__(self, other: Poly) -> Poly:
        result = Poly()
        # We do it both ways to ensure we cover all existing coefficients.
        for x_pow in self.coeff.keys():
            for y_pow, val in self.coeff[x_pow].items():
                result.set_coeff(
                    val+other.get_coeff(x_pow, y_pow), x_pow, y_pow)
        for x_pow in other.coeff.keys():
            for y_pow, val in other.coeff[x_pow].items():
                result.set_coeff(
                    val+self.get_coeff(x_pow, y_pow), x_pow, y_pow)
        return result

    def __sub__(self, other: Poly) -> Poly:
        return self + other * Poly(-1)

    def __mul__(self, other: Poly) -> Poly:
        result = Poly()
        # Cross product. This is slow, especially since we are allocating a new
        # Poly every term but we shouldn't have too many terms.
        for x_pow1 in self.coeff.keys():
            for y_pow1, coeff1 in self.coeff[x_pow1].items():
                for x_pow2 in other.coeff.keys():
                    for y_pow2, coeff2 in other.coeff[x_pow2].items():
                        result += Poly(
                            coeff1*coeff2,
                            x_pow1 + x_pow2,
                            y_pow1+y_pow2,
                        )
        return result

    def __pow__(self, power: int) -> Poly:
        if power < 0:
            raise ArithmeticError("negative powers not supported")
        # We use doubling and multiplication to implement power.
        result = Poly(1)
        doubler = self
        current_pow = 1
        while current_pow <= power:
            if current_pow & power > 0:
                result *= doubler
            doubler *= doubler
            current_pow *= 2
        return result

    def __mod__(self, mod: int) -> Poly:
        result = Poly()
        for x_pow in self.coeff.keys():
            for y_pow, coeff in self.coeff[x_pow].items():
                result.set_coeff(coeff % mod, x_pow, y_pow)
        return result

    def __str__(self) -> str:
        monomials = [
            f"{self._format_num(coeff)} x^{x_pow} y^{y_pow}"
            for x_pow in self.coeff.keys()
            for y_pow, coeff in self.coeff[x_pow].items()
            if coeff != 0
        ]
        if len(monomials) == 0:
            return "0"
        return " + ".join(monomials)

    def __repr__(self) -> str:
        return str(self)

    @staticmethod
    def _format_num(a: int) -> str:
        if -1_000_000 < a < 1_000_000:
            return str(a)
        a_str = str(a)
        sign = ""
        if a_str[0] == "-":
            a_str = a_str[1:]
            sign = "-"
        return f"{sign}{a_str[0]}.{a_str[1:3]}e{len(a_str)-1}"

    def resultant(self, other: Poly) -> Poly:
        self_y = self.highest_y_power()
        other_y = other.highest_y_power()
        sylvester_matrix = PolyMat(self_y + other_y)
        for i in range(other_y):
            for y_pow in range(self_y+1):
                sylvester_matrix[i, i+y_pow] = self._univariate(self_y - y_pow)
        for i in range(self_y):
            for y_pow in range(other_y+1):
                sylvester_matrix[i+other_y, i +
                                 y_pow] = other._univariate(other_y - y_pow)
        return sylvester_matrix.det()

    def _univariate(self, y_pow: int) -> Poly:
        result = Poly()
        for x_pow in range(self.highest_x_power()+1):
            result += Poly(self.get_coeff(x_pow, y_pow), x_pow, 0)
        return result


class PolyMat:
    def __init__(self, dim: int) -> None:
        self._dim = dim
        self._val = [Poly()] * (dim*dim)

    def det(self) -> Poly:
        if self._dim == 1:
            return self._val[0]
        if self._dim == 2:
            return self[0, 0] * self[1, 1] - self[1, 0]*self[0, 1]
        det = Poly()
        for x in range(self._dim):
            det += self[x, 0] * Poly((-1) ** x) * self.submatrix(x, 0).det()
        return det

    def submatrix(self, i: int, j: int) -> PolyMat:
        result = PolyMat(self._dim-1)
        for x in range(self._dim):
            if x == i:
                continue
            sub_x = x if x < i else x-1
            for y in range(self._dim):
                if y == j:
                    continue
                sub_y = y if y < j else y-1
                result[sub_x, sub_y] = self[x, y]
        return result

    @overload
    def __getitem__(self, key: int) -> list[int]:
        pass

    @overload
    def __getitem__(self, key: tuple[int, int]) -> Poly:
        pass

    def __getitem__(self, key: int | tuple[int, int]):
        if isinstance(key, tuple):
            return self._val[key[0]*self._dim + key[1]]
        return self._val[key*self._dim:(key+1)*self._dim]

    def __setitem__(self, key: tuple[int, int], val: Poly) -> None:
        idx = key[0] * self._dim + key[1]
        self._val[idx] = val

    def __str__(self) -> str:
        result = "["
        for i in range(self._dim):
            result += f"\n\t{self[i]}"
        result += "\n]"
        return result


def solve_quadratic(eqn: Poly) -> tuple[float, float]:
    assert eqn.highest_x_power() == 2
    assert eqn.highest_y_power() == 0
    a = eqn.get_coeff(2, 0)
    b = eqn.get_coeff(1, 0)
    c = eqn.get_coeff(0, 0)
    determinant = b*b - 4*a*c
    d = isqrt(determinant)
    if d*d != determinant:
        raise ArithmeticError("Determinant is not perfect square")
    print(a, b, c, d)
    return ((-b + d)/(2*a), (-b-d)/(2*a))


def remove_x_factor(eqn: Poly) -> Poly:
    assert eqn.highest_y_power() == 0
    result = Poly()
    lowest_x: int | None = None
    for x in range(eqn.highest_x_power()+1):
        if eqn.get_coeff(x, 0) == 0:
            continue
        if lowest_x is None:
            lowest_x = x
        result.set_coeff(eqn.get_coeff(x, 0), x-lowest_x, 0)
    return result


# Simple variables for ease of use.
x = Poly(1, 1, 0)
y = Poly(1, 0, 1)
