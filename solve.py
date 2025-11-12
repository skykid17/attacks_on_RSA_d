import decimal
from expression import Poly
from math import gcd, isqrt
from typing import Optional, Sequence
import cmath


def linear(eqn: Poly) -> int:
    assert eqn.highest_x_power() == 1
    assert eqn.highest_y_power() == 0
    assert eqn.highest_z_power() == 0
    c1 = eqn[1, 0, 0]
    c0 = eqn[0, 0, 0]
    if c0 % c1 == 0:
        return -c0 // c1
    raise ArithmeticError("Linear equation has no integer root")


def quadratic(eqn: Poly) -> tuple[int, int]:
    assert eqn.highest_x_power() == 2
    assert eqn.highest_y_power() == 0
    assert eqn.highest_z_power() == 0
    a = eqn[2, 0, 0]
    b = eqn[1, 0, 0]
    c = eqn[0, 0, 0]
    determinant = b*b - 4*a*c
    if determinant < 0:
        raise ArithmeticError("Determinant is negative")
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
    assert eqn.highest_z_power() == 0
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
        highest_coeff = max(x_pow, eqn[x_pow, 0, 0])
    ctx = decimal.Context(prec=(highest_coeff.bit_length()*3 + 20))

    derivative = eqn.derivative_x()
    # Halley's method.
    der_der = derivative.derivative_x()
    # x^(n-1)/x^n is sum of all roots and there are n roots so we use this as
    # our approximation.
    x = (eqn[eqn.highest_x_power() - 1, 0, 0] //
         eqn[eqn.highest_x_power(), 0, 0] //
         eqn.highest_x_power())
    x = x if x > 0 else -x
    if x_guess is not None:
        x = x_guess
    nudge = 0
    while nudge < 1_000:
        f_x = eqn.eval(x, 0, 0)
        if f_x == 0:
            return x
        der_x = ctx.create_decimal(derivative.eval(x, 0, 0))
        der_der_x = ctx.create_decimal(der_der.eval(x, 0, 0))
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
    for x_deviation in range(1_000):
        if eqn.eval(x + x_deviation, 0, 0) == 0:
            return x + x_deviation
        if eqn.eval(x - x_deviation, 0, 0) == 0:
            return x - x_deviation
    raise ArithmeticError("Cannot find actual x.")


def _simplify(eqn: Poly) -> Poly:
    eqn = primitive(eqn)
    # Remove square factors.
    derivative = eqn.derivative_x()
    gcd_with_der = eqn.gcd_x(derivative)
    if gcd_with_der.highest_x_power() > 0:
        return _simplify(eqn / gcd_with_der)
    return eqn


def primitive(eqn: Poly) -> Poly:
    # Simplify polynomial by removing common factors.
    common = 0
    for (coeff, _) in eqn.all_coeffs():
        common = gcd(coeff, common)
    if eqn[eqn.highest_x_power(), 0, 0] < 0:
        # Also flip sign if we are x only.
        common = -common
    return eqn/Poly(common)


def poly_eval(coeffs: list[complex], x: complex) -> complex:
    result = complex(0, 0)
    for c in coeffs:
        result = result * x + c
    return result


def durand_kerner(
        coeffs: Sequence[complex],
        tol: float = 1e-9,
        max_iter: int = 1000,
        verbose: bool = False) -> list[complex]:
    # Convert all coefficients to complex numbers
    c = [complex(c) for c in coeffs]

    # Degree of the polynomial
    n = len(c) - 1
    if n <= 0:
        return []

    # Generate initial root guesses (Aberth-Ehrlich method)
    # Spread the guesses evenly on a circle in the complex plane.
    # The circle's radius is based on the leading coefficients.
    # This is a standard, robust starting point.
    R = (abs(c[-1]) / abs(c[0]))**(1/n)
    v = []
    for k in range(n):
        angle = 2 * cmath.pi * k / n
        # Add a small offset to the angle to avoid symmetries
        v.append(R * cmath.exp(complex(0, angle + 0.4)))

    if verbose:
        print(f"Initial guesses: {v}")

    prev_change = None
    diverge_count = 0
    while True:
        max_change = 0.0
        new_v = list(v)

        for i in range(n):
            # Calculate the update step for root v[i]
            p_vi = poly_eval(c, v[i])
            # Denominator: product(v[i] - v[j] for j != i)
            denominator = complex(1, 0)
            for j in range(n):
                if i == j:
                    continue
                denominator *= (v[i] - v[j])
            if abs(denominator) < 1e-20:
                # This happens with multiple roots or bad guesses.
                # Avoid division by zero and skip update.
                update = complex(0, 0)
            else:
                update = p_vi / denominator
            # Apply the update
            new_v[i] = v[i] - update
            max_change = max(max_change, abs(update))
        v = new_v
        # Check for convergence
        if max_change < tol:
            return v
        # Check for divergence.
        if prev_change is not None and abs(max_change / prev_change) > 1.0:
            diverge_count += 1
            if diverge_count > 1_000:
                print("Roots diverging. Returning known roots.")
                return v
        prev_change = max_change


def find_root_dk(eqn: Poly, x_guess: Optional[int] = None, verbose=False) -> int:
    assert eqn.highest_y_power() == 0
    eqn = _simplify(eqn)
    coeffs = eqn.get_x_coefficients_desc()

    if not coeffs or coeffs[0] == 0:
        raise ArithmeticError("Polynomial is zero.")

    if verbose:
        print(f"Finding roots for degree {len(coeffs)-1} polynomial...")

    # Call durand_kerner
    # This is the line that will likely raise an OverflowError
    potential_roots = durand_kerner(coeffs, tol=1e-9, verbose=verbose)

    if verbose:
        print(f"Got {potential_roots}")

    # Filter for positive integer roots.
    int_roots = []
    for r in potential_roots:
        # Check if imaginary part is tiny
        if abs(r.imag) < 1e-9:
            r_real = r.real
            # Check if it's very close to an integer
            if abs(r_real - round(r_real)) < 1e-9:
                root_val = int(round(r_real))
                if root_val > 0:
                    int_roots.append(root_val)

    int_roots = sorted(list(set(int_roots)))  # Remove duplicates

    if verbose:
        print(f"Found {len(int_roots)
                       } potential positive integer roots: {int_roots}")

    if not int_roots:
        raise ArithmeticError("Durand-Kerner found no positive integer roots.")

    # Return the "best" root.
    # pick the one closest to x_guess, or just the smallest.
    if x_guess is not None:
        best_root = min(int_roots, key=lambda r: abs(r - x_guess))
        return best_root

    return int_roots[0]  # Return the smallest positive integer root
