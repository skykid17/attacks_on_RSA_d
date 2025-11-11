import decimal
from expression import Poly
from math import gcd, isqrt
from typing import Optional


def linear(eqn: Poly) -> int:
    assert eqn.highest_x_power() == 1
    assert eqn.highest_y_power() == 0
    c1 = eqn.get_coeff(1, 0)
    c0 = eqn.get_coeff(0, 0)
    if c0 % c1 == 0:
        return -c0 // c1
    raise ArithmeticError("Linear equation has no integer root")


def quadratic(eqn: Poly) -> tuple[int, int]:
    assert eqn.highest_x_power() == 2
    assert eqn.highest_y_power() == 0
    a = eqn.get_coeff(2, 0)
    b = eqn.get_coeff(1, 0)
    c = eqn.get_coeff(0, 0)
    determinant = b*b - 4*a*c
    d = isqrt(determinant)
    if d*d != determinant:
        raise ArithmeticError("Determinant is not perfect square")
    if (d-b) % (2*a) != 0:
        raise ArithmeticError("-b+d is not a multiple of 2*a")
    if (-b-d) % (2*a) != 0:
        raise ArithmeticError("-b-d is not a multiple of 2*a")
    return ((-b + d)//(2*a), (-b-d)//(2*a))


def first_root(eqn: Poly, verbose=False, x_guess: Optional[int] = None) -> int:
    assert eqn.highest_y_power() == 0
    assert eqn.highest_x_power() > 0

    eqn = _simplify(eqn)
    if verbose:
        print(f"After simplification: {eqn}")

    if eqn.highest_x_power() == 1:
        return linear(eqn)
    if eqn.highest_x_power() == 2:
        x1, x2 = quadratic(eqn)
        x1, x2 = min(x1, x2), max(x1, x2)
        if x1 > 0:
            return x1
        return x2

    highest_coeff = 0
    for x_pow in range(eqn.highest_x_power()+1):
        highest_coeff = max(x_pow, eqn.get_coeff(x_pow, 0))
    ctx = decimal.Context(prec=(highest_coeff.bit_length()*3 + 20))

    derivative = eqn.derivative_x()
    # Halley's method.
    der_der = derivative.derivative_x()
    # x^(n-1)/x^n is sum of all roots and there are n roots so we use this as
    # our approximation.
    x = (eqn.get_coeff(eqn.highest_x_power() - 1, 0) //
         eqn.get_coeff(eqn.highest_x_power(), 0) //
         eqn.highest_x_power())
    x = x if x > 0 else -x
    if x_guess is not None:
        x = x_guess
    nudge = 0
    while nudge < 1_000:
        f_x = eqn.eval(x, 0)
        if f_x == 0:
            return x
        der_x = ctx.create_decimal(derivative.eval(x, 0))
        der_der_x = ctx.create_decimal(der_der.eval(x, 0))
        denominator = der_x*der_x - der_x * der_der_x / 2
        if denominator == 0:
            # We stop if we're stuck nudging.
            nudge += 1
            x += 1
            continue
        update_amount = int(((f_x * der_x) / denominator).to_integral_value())
        if update_amount == 0:
            # We have hit the best we can do. We should be in the
            # neighbourhood.
            break
        x = x-update_amount
    for x_deviation in range(10_000):
        if eqn.eval(x + x_deviation, 0) == 0:
            return x + x_deviation
        if eqn.eval(x - x_deviation, 0) == 0:
            return x - x_deviation
    raise ArithmeticError("Cannot find actual x.")


def _simplify(eqn: Poly) -> Poly:
    # Simplify polynomial by removing common factors.
    common = 0
    for x_pow in range(eqn.highest_x_power()+1):
        common = gcd(abs(eqn.get_coeff(x_pow, 0)), common)
    if eqn.get_coeff(eqn.highest_x_power(), 0) < 0:
        common = -common
    new_eqn = Poly()
    for x_pow in range(eqn.highest_x_power()+1):
        new_eqn.set_coeff(eqn.get_coeff(x_pow, 0) // common, x_pow, 0)

    # Remove square factors.
    derivative = new_eqn.derivative_x()
    gcd_with_der = new_eqn.gcd_x(derivative)
    if gcd_with_der.highest_x_power() > 0:
        return _simplify(new_eqn / gcd_with_der)
    return new_eqn
