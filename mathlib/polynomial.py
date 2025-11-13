from __future__ import annotations
from collections import defaultdict
from math import gcd
import sys
from typing import Generator, overload
type Pow = tuple[int, int, int]


class Poly:
    """
    Poly is a simple polynomial that stores a sum of a*x^b*y^c terms.
    """

    FULL_THRESHOLD = 10 ** 6

    coeff: defaultdict[int, defaultdict[int, defaultdict[int, int]]]

    def __init__(self, constant: int = 0, pow: Pow = (0, 0, 0)) -> None:
        self.coeff = defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(int)
            ))
        self[pow] = constant

    def is_zero(self) -> bool:
        return len(self.coeff.keys()) == 0

    def all_coeffs(self) -> Generator[tuple[int, Pow]]:
        for x_pow, coeff_x in self.coeff.items():
            for y_pow, coeff_y in coeff_x.items():
                for z_pow, coeff in coeff_y.items():
                    if coeff != 0:
                        yield (coeff, (x_pow, y_pow, z_pow))

    def __getitem__(self, pow: Pow) -> int:
        x_pow, y_pow, z_pow = pow
        return self.coeff[x_pow][y_pow][z_pow]

    def __setitem__(self, pow: Pow, coeff: int) -> None:
        x_pow, y_pow, z_pow = pow
        if coeff == 0:
            if (x_pow not in self.coeff or
                y_pow not in self.coeff[x_pow] or
                    z_pow not in self.coeff[x_pow][y_pow]):
                return
            # Delete to save a bit of memory.
            del self.coeff[x_pow][y_pow][z_pow]
            if len(self.coeff[x_pow][y_pow]) == 0:
                del self.coeff[x_pow][y_pow]
            if len(self.coeff[x_pow]) == 0:
                del self.coeff[x_pow]
            return
        self.coeff[x_pow][y_pow][z_pow] = coeff

    def highest_x_power(self) -> int:
        return max(self.coeff.keys(), default=0)

    def highest_y_power(self) -> int:
        return max(
            (y_pow
             for x_pow in self.coeff.keys()
             for y_pow in self.coeff[x_pow].keys()),
            default=0,
        )

    def highest_z_power(self) -> int:
        return max(
            (z_pow
             for (_, (_, _, z_pow))
             in self.all_coeffs()),
            default=0
        )

    def derivative_x(self) -> Poly:
        result = Poly()
        for x_pow in self.coeff.keys():
            if x_pow == 0:
                # Constant is derived away.
                continue
            for y_pow, y_val in self.coeff[x_pow].items():
                for z_pow, val in y_val.items():
                    result[x_pow-1, y_pow, z_pow] += val * x_pow
        return result

    def eval(self, x: int, y: int, z: int) -> int:
        result = 0
        x_val = 1
        for x_pow in range(self.highest_x_power()):
            for y_pow, y_val in self.coeff[x_pow].items():
                y_mul = y ** y_pow
                for z_pow, val in y_val.items():
                    result += x_val * y_mul * (z ** z_pow) * val
            x_val *= x
        # For the last cycle, don't bother calculating x ** (highest+1).
        for y_pow, y_val in self.coeff[self.highest_x_power()].items():
            y_mul = y ** y_pow
            for z_pow, val in y_val.items():
                result += x_val * y_mul * (z ** z_pow) * val
        return result

    def eval_poly(self, x: Poly | int, y: Poly | int, z: Poly | int) -> Poly:
        result = Poly()
        for x_pow in self.coeff.keys():
            x_mul = x ** x_pow
            for y_pow, y_val in self.coeff[x_pow].items():
                y_mul = y ** y_pow
                for z_pow, val in y_val.items():
                    z_mul = z ** z_pow
                    result += x_mul * y_mul * z_mul * val
        return result

    def sub_xy(self, xy: Poly) -> Poly:
        result = Poly()
        for x_pow in self.coeff.keys():
            for y_pow, y_val in self.coeff[x_pow].items():
                xy_replaced = min(x_pow, y_pow)
                xy_mul = xy ** xy_replaced
                remaining_mul = Poly(
                    1, (x_pow-xy_replaced, y_pow-xy_replaced, 0))
                new_xy_mul = xy_mul * remaining_mul
                for z_pow, val in y_val.items():
                    z_mul = Poly(1, (0, 0, z_pow))
                    result += new_xy_mul * z_mul * Poly(val)
        return result

    def copy(self) -> Poly:
        result = Poly()
        for (coeff, pow) in self.all_coeffs():
            result[pow] = coeff
        return result

    def __add__(self, other: Poly | int) -> Poly:
        result = self.copy()
        if isinstance(other, int):
            result[0, 0, 0] += other
            return result
        for (coeff, pow) in other.all_coeffs():
            result[pow] += coeff
        return result

    def __sub__(self, other: Poly | int) -> Poly:
        if isinstance(other, Poly):
            return self + other * Poly(-1)
        result = self.copy()
        result[0, 0, 0] -= other
        return result

    def __neg__(self) -> Poly:
        return self * -1

    def __rsub__(self, other: int) -> Poly:
        return self * -1 + other

    def __mul__(self, other: Poly | int) -> Poly:
        result = Poly()
        if isinstance(other, int):
            for (coeff, pow) in self.all_coeffs():
                result[pow] = coeff * other
            return result
        # Cross product. This is slow, especially since we are allocating a new
        # Poly every term but we shouldn't have too many terms.
        for (coeff1, (x_pow1, y_pow1, z_pow1)) in self.all_coeffs():
            for (coeff2, (x_pow2, y_pow2, z_pow2)) in other.all_coeffs():
                new_pow = (
                    x_pow1 + x_pow2,
                    y_pow1 + y_pow2,
                    z_pow1 + z_pow2,
                )
                result[new_pow] += coeff1*coeff2
        return result

    def __rmul__(self, other: int) -> Poly:
        return self * other

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
            if current_pow == power:
                break
            doubler *= doubler
            current_pow *= 2
        return result

    def __mod__(self, mod: int) -> Poly:
        result = Poly()
        for (coeff, pow) in self.all_coeffs():
            result[pow] = coeff % mod
        return result

    def __truediv__(self, divisor: Poly | int) -> Poly:
        assert self.highest_z_power() == 0
        if isinstance(divisor, int):
            divisor = Poly(divisor, (0, 0, 0))
        assert divisor.highest_z_power() == 0

        result = Poly()
        divisor_y_pow = divisor.highest_y_power()
        highest_y = divisor.univariate_y(divisor_y_pow)
        numerator = self
        while not numerator.is_zero():
            y_pow = numerator.highest_y_power()
            if y_pow < divisor_y_pow:
                raise ArithmeticError(
                    f"cannot divide {self} exactly by {divisor}: " +
                    f"remainder {numerator}")
            num_y = numerator.univariate_y(y_pow)
            k, divided, rem = num_y._div_x(highest_y)
            if k > 1 or not rem.is_zero():
                raise ArithmeticError(
                    f"cannot divide {self} exactly by {divisor}: " +
                    f"remainder {numerator}")
            factor = divided * Poly(1, (0, y_pow - divisor_y_pow, 0))
            result += factor
            numerator -= factor*divisor
        return result

    def gcd_x(self, other: Poly) -> Poly:
        assert self.highest_y_power() == 0
        assert other.highest_y_power() == 0
        a_content, a = primitive_content(self)
        b_content, b = primitive_content(other)
        while not a.is_zero():
            # b._div_x(a) returns the result of divmod(kb, a).
            # This means that if k > 1, we are actually calculating
            # gcd(ka, b) where k is not a factor of a.
            # As such, we can pull out all content from remainder as it was
            # primitive all along.
            _, _, rem = b._div_x(a, rem_only=True)
            a, b = primitive(rem), a
        return gcd(a_content, b_content) * b

    def _div_x(self,
               divisor: Poly,
               rem_only: bool = False) -> tuple[int, Poly, Poly]:
        """
        _div_x divides two univariate X polynomials. If required, it will
        multiply the numerator by a constant, this is returned as an integer.
        The second and third values are the result and remainder of the
        division.
        """
        assert self.highest_y_power() == 0
        assert self.highest_z_power() == 0
        assert divisor.highest_y_power() == 0
        assert divisor.highest_z_power() == 0
        numerator = self
        result = Poly()
        lead_den = divisor[divisor.highest_x_power(), 0, 0]
        dividend_mul = 1
        assert lead_den != 0
        while not numerator.is_zero():
            num_x_pow = numerator.highest_x_power()
            div_x_pow = divisor.highest_x_power()
            if num_x_pow < div_x_pow:
                # Not fully divisible, return.
                return dividend_mul, result, numerator
            lead_num = numerator[numerator.highest_x_power(), 0, 0]
            quo, rem = divmod(lead_num, lead_den)
            if rem != 0:
                # Not fully divisible. We multiply by a factor to continue.
                mul_factor = lead_den // gcd(rem, lead_den)
                numerator *= mul_factor
                dividend_mul *= mul_factor
                if not rem_only:
                    result *= mul_factor
                continue
            factor = Poly(quo, (num_x_pow-div_x_pow, 0, 0))
            numerator -= factor * divisor
            if not rem_only:
                result += factor
        return dividend_mul, result, numerator

    def __str__(self) -> str:
        monomials = [
            self._format_mono(coeff, x_pow, y_pow, z_pow)
            for x_pow, x_val in self.coeff.items()
            for y_pow, y_val in x_val.items()
            for z_pow, coeff in y_val.items()
            if coeff != 0
        ]
        if len(monomials) == 0:
            return "0"
        return " + ".join(monomials)

    def __repr__(self) -> str:
        return str(self)

    @staticmethod
    def _format_mono(a: int, x_pow: int, y_pow: int, z_pow: int) -> str:
        if x_pow == 0:
            x_suffix = ""
        elif x_pow == 1:
            x_suffix = "x"
        else:
            x_suffix = f"x^{x_pow} "
        if y_pow == 0:
            y_suffix = ""
        elif y_pow == 1:
            y_suffix = "y"
        else:
            y_suffix = f"y^{y_pow} "
        if z_pow == 0:
            z_suffix = ""
        elif z_pow == 1:
            z_suffix = "z"
        else:
            z_suffix = f"z^{z_pow} "
        return f"{Poly._format_num(a)}{x_suffix}{y_suffix}{z_suffix}"

    @staticmethod
    def _format_num(a: int) -> str:
        # We need to override the limits to actually print it out...
        original_limit = sys.get_int_max_str_digits()
        sys.set_int_max_str_digits(1_000_000)
        a_str = str(a)
        sys.set_int_max_str_digits(original_limit)
        if -Poly.FULL_THRESHOLD < a < Poly.FULL_THRESHOLD:
            return a_str
        sign = ""
        if a_str[0] == "-":
            a_str = a_str[1:]
            sign = "-"
        return f"{sign}{a_str[0]}.{a_str[1:3]}e{len(a_str)-1}"

    def univariate_y(self, y_pow: int) -> Poly:
        result = Poly()
        for (coeff, (x_pow, y, z_pow)) in self.all_coeffs():
            if y == y_pow:
                result[x_pow, 0, z_pow] = coeff
        return result

    def has_y(self) -> bool:
        """Checks if the polynomial has any y-terms."""
        return self.highest_y_power() > 0

    def get_x_coefficients_desc(self) -> list[int]:
        """
        Gets coefficients for a univariate-in-x polynomial in ascending order
        [c0, c1, c2, ..., cn]
        """
        assert self.highest_y_power() == 0, "Not a univariate-in-x polynomial"
        assert self.highest_z_power() == 0, "Not a univariate-in-x polynomial"

        high_x = self.highest_x_power()
        return [self[i, 0, 0] for i in range(high_x, -1, -1)]

    def degree(self) -> int:
        """Returns the highest total degree of any term (x_pow + y_pow)."""
        max_deg = 0
        for x_pow, y_map in self.coeff.items():
            for y_pow in y_map.keys():
                max_deg = max(max_deg, x_pow + y_pow)
        return max_deg

def primitive_content(eqn: Poly) -> tuple[int, Poly]:
    if eqn.is_zero():
        return 0, eqn
    # Simplify polynomial by removing common factors.
    common = 0
    for (coeff, _) in eqn.all_coeffs():
        common = gcd(coeff, common)
    if eqn[eqn.highest_x_power(), 0, 0] < 0:
        # Also flip sign if we are x only.
        common = -common
    return common, eqn/common


def primitive(eqn: Poly) -> Poly:
    _, p = primitive_content(eqn)
    return p



class PolyMat:
    def __init__(self, dim: int) -> None:
        assert dim > 0
        self._dim = dim
        self._val = [Poly()] * (dim*dim)

    def det(self) -> Poly:
        if self._dim == 1:
            return self._val[0]
        if self._dim == 2:
            return self[0, 0] * self[1, 1] - self[1, 0]*self[0, 1]
        # Bareiss algorithm.
        mat = PolyMat(self._dim)  # Make a copy.
        for i in range(self._dim):
            mat[i] = self[i]
        prev_pivot = Poly(1)
        sign = 1
        for k in range(0, self._dim-1):
            pivot_row = k
            while mat[pivot_row, k].is_zero():
                pivot_row += 1
                if pivot_row == self._dim:
                    # No non-zero pivot, the determinant is zero.
                    return Poly(0)
            if pivot_row != k:
                # Need to do the swap.
                # Switching rows force us to switch signs.
                mat[pivot_row], mat[k] = mat[k], mat[pivot_row]
                sign = -sign
            pivot = mat[k, k]
            for i in range(k+1, self._dim):
                for j in range(k+1, self._dim):
                    num = mat[k, k]*mat[i, j] - mat[i, k]*mat[k, j]
                    mat[i, j] = num / prev_pivot
                mat[i, k] = Poly(0)
            prev_pivot = pivot
        return Poly(sign) * mat[self._dim-1, self._dim-1]

    @overload
    def __getitem__(self, key: int) -> list[Poly]:
        pass

    @overload
    def __getitem__(self, key: tuple[int, int]) -> Poly:
        pass

    def __getitem__(self, key: int | tuple[int, int]):
        if isinstance(key, tuple):
            assert key[0] < self._dim
            assert key[1] < self._dim
            return self._val[key[0]*self._dim + key[1]]
        assert key < self._dim
        return self._val[key*self._dim:(key+1)*self._dim]

    @overload
    def __setitem__(self, key: int, val: list[Poly]) -> None:
        pass

    @overload
    def __setitem__(self, key: tuple[int, int], val: Poly) -> None:
        pass

    def __setitem__(self,
                    key: int | tuple[int, int],
                    val: list[Poly] | Poly) -> None:
        if isinstance(key, tuple) and isinstance(val, Poly):
            assert key[0] >= 0 and key[0] < self._dim
            assert key[1] >= 0 and key[1] < self._dim
            idx = key[0] * self._dim + key[1]
            self._val[idx] = val
            return
        if isinstance(key, int) and isinstance(val, list):
            assert key >= 0 and key < self._dim
            assert len(val) == self._dim
            for i, v in enumerate(val):
                self[key, i] = v
            return
        raise AssertionError("invalid types")

    def __str__(self) -> str:
        result = "["
        for i in range(self._dim):
            result += f"\n\t{self[i]}"
        result += "\n]"
        return result
