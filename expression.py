from __future__ import annotations
from collections import defaultdict
import sys
from typing import overload
from fractions import Fraction


class Poly:
    """
    Poly is a simple polynomial that stores a sum of a*x^b*y^c terms.
    """

    FULL_THRESHOLD = 10 ** 6

    def __init__(
            self,
            constant: int = 0, x_pow: int = 0, y_pow: int = 0) -> None:
        self.coeff: defaultdict[int, defaultdict[int, int]] = defaultdict(
            lambda: defaultdict(int))
        self.set_coeff(constant, x_pow, y_pow)

    def is_zero(self) -> bool:
        return len(self.coeff.keys()) == 0

    def all_coeffs(self) -> list[tuple[int, int, int]]:
        return [
            (coeff, x, y)
            for x in self.coeff.keys()
            for y, coeff in self.coeff[x].items()
            if coeff != 0
        ]

    def get_coeff(self, x_pow: int, y_pow: int) -> int:
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
        return max(self.coeff.keys(), default=0)

    def highest_y_power(self) -> int:
        return max(
            (y_pow for x_pow in self.coeff.keys()
             for y_pow in self.coeff[x_pow].keys()),
            default=0,
        )

    def derivative_x(self) -> Poly:
        result = Poly()
        for x_pow in self.coeff.keys():
            if x_pow == 0:
                # Constant is derived away.
                continue
            multiplier = Poly(x_pow) * (x ** (x_pow-1))
            for y_pow, val in self.coeff[x_pow].items():
                result += multiplier * (y ** y_pow) * Poly(val)
        return result

    def eval(self, x: int, y: int) -> int:
        result = 0
        x_val = 1
        for x_pow in range(self.highest_x_power()):
            for y_pow, val in self.coeff[x_pow].items():
                result += x_val * (y ** y_pow) * val
            x_val *= x
        # For the last cycle, don't bother calculating x ** (highest+1).
        for y_pow, val in self.coeff[self.highest_x_power()].items():
            result += x_val * (y ** y_pow) * val
        return result

    def eval_poly(self, x: Poly, y: Poly) -> Poly:
        result = Poly()
        for x_pow in self.coeff.keys():
            for y_pow, val in self.coeff[x_pow].items():
                result += (x ** x_pow) * (y ** y_pow) * Poly(val)
        return result

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
                            x_pow1+x_pow2,
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
            if current_pow == power:
                break
            doubler *= doubler
            current_pow *= 2
        return result

    def __mod__(self, mod: int) -> Poly:
        result = Poly()
        for x_pow in self.coeff.keys():
            for y_pow, coeff in self.coeff[x_pow].items():
                result.set_coeff(coeff % mod, x_pow, y_pow)
        return result

    def __truediv__(self, divisor: Poly) -> Poly:
        result = Poly()
        divisor_y_pow = divisor.highest_y_power()
        highest_y = divisor._univariate(divisor_y_pow)
        numerator = self
        while not numerator.is_zero():
            y_pow = numerator.highest_y_power()
            if y_pow < divisor_y_pow:
                raise ArithmeticError(
                    f"cannot divide {self} exactly by {divisor}: " +
                    f"remainder {numerator}")
            num_y = numerator._univariate(y_pow)
            divided, rem = num_y._div_x(highest_y)
            if not rem.is_zero():
                raise ArithmeticError(
                    f"cannot divide {self} exactly by {divisor}: " +
                    f"remainder {numerator}")
            factor = divided * (y ** (y_pow - divisor_y_pow))
            result += factor
            numerator -= factor*divisor
        return result

    def gcd_x(self, other: Poly) -> Poly:
        assert self.highest_y_power() == 0
        assert other.highest_y_power() == 0
        a, b = self, other
        while not a.is_zero():
            _, rem = b._div_x(a)
            if rem.highest_x_power() >= a.highest_x_power():
                # They do not share an integer gcd.
                return Poly(1)
            a, b = rem, a
        return b

    def _div_x(self, divisor: Poly) -> tuple[Poly, Poly]:
        assert self.highest_y_power() == 0
        assert divisor.highest_y_power() == 0
        numerator = self
        result = Poly()
        lead_den = divisor.get_coeff(divisor.highest_x_power(), 0)
        assert lead_den != 0
        while not numerator.is_zero():
            num_x_pow = numerator.highest_x_power()
            div_x_pow = divisor.highest_x_power()
            if num_x_pow < div_x_pow:
                # Not fully divisible, return.
                return result, numerator
            lead_num = numerator.get_coeff(numerator.highest_x_power(), 0)
            if abs(lead_num) < abs(lead_den):
                # Not fully divisible, return.
                return result, numerator
            factor = Poly(lead_num // lead_den) * (x ** (num_x_pow-div_x_pow))
            numerator -= factor * divisor
            result += factor
        assert numerator.is_zero()
        return result, numerator

    def __str__(self) -> str:
        monomials = [
            self._format_mono(coeff, x_pow, y_pow)
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
    def _format_mono(a: int, x_pow: int, y_pow: int) -> str:
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
        return f"{Poly._format_num(a)}{x_suffix}{y_suffix}"

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
    
    # --- Add this method INSIDE the Poly class in expression.py ---

    def get_x_coefficients_desc(self) -> list[int]:
        """
        Gets coefficients for a univariate-in-x polynomial in descending order
        [c_n, c_{n-1}, ..., c_0]
        """
        assert self.highest_y_power() == 0, "Not a univariate-in-x polynomial"
        
        high_x = self.highest_x_power()
        # Create a list of zeros [c_0, c_1, ..., c_n]
        coeffs_asc = [0] * (high_x + 1)
        
        for x_pow, y_map in self.coeff.items():
            if 0 in y_map: # Check if y_pow=0 exists
                coeffs_asc[x_pow] = y_map[0]
        
        # Reverse to get [c_n, c_{n-1}, ..., c_0]
        return list(reversed(coeffs_asc))
    def has_y(self) -> bool:
        """Checks if the polynomial has any y-terms."""
        return self.highest_y_power() > 0

    def get_x_coefficients_asc(self) -> list[int]:
        """
        Gets coefficients for a univariate-in-x polynomial in ascending order
        [c0, c1, c2, ..., cn]
        """
        assert self.highest_y_power() == 0, "Not a univariate-in-x polynomial"
        
        high_x = self.highest_x_power()
        # Create a list of zeros [c0, c1, ..., cn]
        coeffs_asc = [0] * (high_x + 1)
        
        for x_pow, y_map in self.coeff.items():
            if 0 in y_map: # Check if y_pow=0 exists
                coeffs_asc[x_pow] = y_map[0]
        
        return coeffs_asc

    def degree(self) -> int:
        """Returns the highest total degree of any term (x_pow + y_pow)."""
        max_deg = 0
        for x_pow, y_map in self.coeff.items():
            for y_pow in y_map.keys():
                max_deg = max(max_deg, x_pow + y_pow)
        return max_deg

class PolyMat:
    def __init__(self, dim: int) -> None:
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
            return self._val[key[0]*self._dim + key[1]]
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
