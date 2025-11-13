from .complex import Complex, ComplexContext
from .polynomial import Poly, PolyMat, primitive, primitive_content
from decimal import MAX_EMAX
from math import gcd, isqrt
from typing import Optional, Sequence


# Simple variables for ease of use.
x = Poly(1, (1, 0, 0))
y = Poly(1, (0, 1, 0))
z = Poly(1, (0, 0, 1))


def remove_x_factor(eqn: Poly) -> Poly:
    assert eqn.highest_y_power() == 0
    assert eqn.highest_z_power() == 0
    result = Poly()
    lowest_x: int | None = None
    for x in range(eqn.highest_x_power()+1):
        coeff = eqn[x, 0, 0]
        if coeff == 0:
            continue
        if lowest_x is None:
            lowest_x = x
        result[x-lowest_x, 0, 0] = coeff
    return result


def resultant(a: Poly, b: Poly) -> Poly:
    self_y = a.highest_y_power()
    other_y = b.highest_y_power()
    if self_y + other_y == 0:
        return Poly(0)
    sylvester_matrix = PolyMat(self_y + other_y)
    for i in range(other_y):
        for y_pow in range(self_y+1):
            sylvester_matrix[i, i+y_pow] = a.univariate_y(self_y - y_pow)
    for i in range(self_y):
        for y_pow in range(other_y+1):
            sylvester_matrix[i+other_y, i +
                             y_pow] = b.univariate_y(other_y - y_pow)
    return sylvester_matrix.det()


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
    return ((-b - d)//(2*a), (-b+d)//(2*a))




def _simplify(eqn: Poly) -> Poly:
    eqn = primitive(eqn)
    # Remove square factors.
    derivative = eqn.derivative_x()
    gcd_with_der = eqn.gcd_x(derivative)
    if gcd_with_der.highest_x_power() > 0:
        return _simplify(eqn / gcd_with_der)
    return eqn


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
    root_scale = ctx.create(coeffs[-1] // coeffs[0]) ** (1/n)
    v = [ctx.create(0.4, 0.9) ** i * root_scale for i in range(n)]

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
    x_pow = eqn.highest_x_power()
    if x_pow == 1:
        return linear(eqn)
    if x_pow == 2:
        quad_roots = quadratic(eqn)
        if x_guess is not None:
            best_root = min(quad_roots, key=lambda r: abs(r - x_guess))
            return best_root
        return quad_roots[0]  # Return the smaller root.

    eqn = _simplify(eqn)
    coeffs = eqn.get_x_coefficients_desc()
    der_coeffs = eqn.derivative_x().get_x_coefficients_desc()

    if not coeffs or coeffs[0] == 0:
        raise ArithmeticError("Polynomial is zero.")

    root_scale = (len(str(coeffs[-1])) - len(str(coeffs[0]))) // len(coeffs)
    prec = root_scale * 2 + 50

    if verbose:
        print(f"Finding roots for degree {len(coeffs)-1} polynomial...")
        print(f"Using {prec=}")

    ctx = ComplexContext(
        prec=prec,
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
            print(f"Found {len(int_roots)} potential positive integer roots: {int_roots}")

        # Return the "best" root.
        # pick the one closest to x_guess, or just the smallest.
        if x_guess is not None:
            best_root = min(int_roots, key=lambda r: abs(r - x_guess))
            return best_root

        return int_roots[0]  # Return the smallest positive integer root
