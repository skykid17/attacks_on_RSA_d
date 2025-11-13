from decimal import MAX_EMAX
from expression import Complex, ComplexContext, Poly
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
    return eqn/common


def aberth_ehrlich(
        ctx: ComplexContext,
        coeffs: Sequence[int],
        der_coeffs: Sequence[int],
        tol: float = 1e-9) -> list[Complex]:
    # Convert all coefficients to complex numbers
    c = [ctx.create(c)/ctx.create(coeffs[0]) for c in coeffs]
    der_c = [ctx.create(c)/ctx.create(coeffs[0]) for c in der_coeffs]
    one = ctx.create(1)

    # Degree of the polynomial
    n = len(c) - 1
    if n <= 0:
        return []

    # Generate initial root guesses. Spread the guesses out.
    v = [ctx.create(0.4, 0.9) ** i for i in range(n)]

    def poly_eval(coeffs: Sequence[Complex],
                  x: Complex) -> Complex:
        result = ctx.create(0)
        for c in coeffs:
            result = result * x + c
        return result

    prev_change = None
    diverge_count = 0
    while True:
        max_change = ctx.create_decimal(0)

        for i in range(n):
            # Calculate the update step for root v[i]
            p_vi = poly_eval(c, v[i])
            p_der_vi = poly_eval(der_c, v[i])
            denominator = ctx.create(0)
            for j in range(n):
                if i == j:
                    continue
                denominator += one / (v[i] - v[j])
            update = p_vi / (p_der_vi - p_vi * denominator)
            v[i] -= update
            max_change = max(max_change, abs(update))
        # Check for convergence
        if max_change < tol:
            return v
        # Check for divergence.
        if (prev_change is not None and
                abs(max_change) / abs(prev_change) >= 1.0):
            diverge_count += 1
            if diverge_count > 10_000:
                print("Roots diverging. Giving up.")
                return []
        prev_change = max_change


def first_root(eqn: Poly, x_guess: Optional[int] = None, verbose=False) -> int:
    assert eqn.highest_y_power() == 0
    eqn = _simplify(eqn)
    coeffs = eqn.get_x_coefficients_desc()
    der_coeffs = eqn.derivative_x().get_x_coefficients_desc()

    if not coeffs or coeffs[0] == 0:
        raise ArithmeticError("Polynomial is zero.")

    if verbose:
        print(f"Finding roots for degree {len(coeffs)-1} polynomial...")

    root_scale = (len(str(coeffs[-1])) - len(str(coeffs[0]))) // len(coeffs)
    ctx = ComplexContext(
        prec=root_scale * 2 + 50,
        Emax=MAX_EMAX,
    )
    with ctx:
        potential_roots = aberth_ehrlich(ctx, coeffs, der_coeffs, tol=1e-9)

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
                        assert eqn.eval(root_val, 0, 0) == 0
                        int_roots.append(root_val)

        int_roots = sorted(list(set(int_roots)))  # Remove duplicates

        if verbose:
            print(f"Found {len(int_roots)
                           } potential positive integer roots: {int_roots}")

        if not int_roots:
            raise ArithmeticError(
                "Durand-Kerner found no positive integer roots.")

        # Return the "best" root.
        # pick the one closest to x_guess, or just the smallest.
        if x_guess is not None:
            best_root = min(int_roots, key=lambda r: abs(r - x_guess))
            return best_root

        return int_roots[0]  # Return the smallest positive integer root
