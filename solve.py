from expression import Poly
from math import isqrt


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


def first_root(eqn: Poly) -> int:
    assert eqn.highest_y_power() == 0
    assert eqn.highest_x_power() > 0

    # Newton method.
    derivative = eqn.derivative_x()
    x = 100
    if eqn.get_coeff(1, 0) != 0:
        # Initial best guess based on -A0/A1.
        x = max(x, -eqn.get_coeff(0, 0) // eqn.get_coeff(1, 0))
    while True:
        f_x = eqn.eval(x, 0)
        der_x = derivative.eval(x, 0)
        if der_x == 0:
            raise ArithmeticError("Derivative hit zero.")
        update_amount = f_x // der_x
        if update_amount == 0:
            # We have hit the best we can do. We should be in the
            # neighbourhood.
            break
        x = x-update_amount
        if x < 0:
            raise ArithmeticError("x going negative.")
    for x_deviation in range(1_000_000):
        if eqn.eval(x + x_deviation, 0) == 0:
            return x + x_deviation
        if eqn.eval(x - x_deviation, 0) == 0:
            return x - x_deviation
    raise ArithmeticError("Cannot find actual x.")
