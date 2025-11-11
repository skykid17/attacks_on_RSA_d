from expression import Poly, x, y, remove_x_factor
from math import ceil
from typing import Optional, Tuple
from solve import first_root, quadratic
from fpylll import LLL, IntegerMatrix


def boneh_durfee_attack(
    n: int,
    e: int,
    m: int = 7,
    t: int = 3,
    X_bound: Optional[int] = None,
    Y_bound: Optional[int] = None,
    delta: float = 0.18,
    timeout_sec: float = 60,
    verbose=False,
) -> Optional[Tuple[int, int, int]]:
    if verbose:
        print(f"Solving n={Poly._format_num(n)} e={Poly._format_num(e)}")
    if n % 2 == 0:
        # n must be even, as such, we can already factorize it!
        phi = n // 2 - 1
        return n // 2, 2, pow(e, -1, phi)

    # Some default values. Y should be in the range p+q and X should be in the
    # range e^delta.
    if X_bound is None:
        X_bound = 2 ** int(ceil(n.bit_length() * delta))
    if Y_bound is None:
        Y_bound = 2 ** (n.bit_length() // 2 + 1)

    # We are trying to solve:
    # f(x, y) = x(A-y)+1 = 0 (mod e).
    # A = N+1, y = p+q
    #
    # This is derived from:
    # ed = k phi(n) + 1
    # ed = k(N-p-q+1) + 1
    # k(N-p-q+1) + 1 = 0 (mod e).
    A = n + 1
    f = x * (Poly(A) - y) + Poly(1)
    e_poly = Poly(e)

    if verbose:
        print(f"Function we are trying to solve is (mod e={
              Poly._format_num(e)}):")
        print(f)
        print()

    polynomials: list[Poly] = []
    # Boneh and Durfee found that the best m is 7 and t is 3.
    # That means our equations are 0 (mod e^7) and have a maximum power of y^3.
    f_pow: list[Poly] = [Poly(1), f]
    fp = f
    for _ in range(2, m+1):
        fp *= f
        f_pow.append(fp)
    for k in range(m + 1):
        # x-shifts: x^i * f(x, y)^k * e^(m-k) = 0 (mod e^m)
        # For each k, we use g_ik for i = 0, ..., m-k.
        for i in range(m-k+1):
            g_ik = (x**i) * f_pow[k] * (e_poly ** (m-k))
            polynomials.append(g_ik)
    for k in range(m + 1):
        # y-shifts: y^i * f(x, y)^k * e^(m-k) = 0 (mod e^m)
        # For each k, we use h_lk for l = 1, ..., t.
        for l_pow in range(1, t+1):
            h_lk = (y**l_pow) * f_pow[k] * (e_poly ** (m-k))
            polynomials.append(h_lk)

    if verbose:
        print("Equations before LLL:")
        for i in range(4):
            print(polynomials[i])
        print()
        print("Doing LLL...")
        print()
    polynomials = _do_lll(polynomials, X_bound, Y_bound)
    if verbose:
        print("Equations after LLL:")
        for i in range(4):
            print(polynomials[i])
        print()
    root_x = None
    for i in range(2, len(polynomials)):
        poly_i = polynomials[i]
        for j in range(1, i):
            poly_j = polynomials[j]
            r_candidate = poly_i.resultant(poly_j)
            r_candidate = remove_x_factor(r_candidate)
            if r_candidate.highest_x_power() == 0:
                # Useless resultant. Skip.
                continue
            if verbose:
                print(f"Using {i=} {j=}")
                print("Chosen resultant:")
                print(r_candidate)
            try:
                eqn_1 = poly_i
                eqn_2 = poly_j
                root_x = first_root(
                    remove_x_factor(r_candidate),
                    x_guess=X_bound//2,
                    verbose=verbose)
                if verbose:
                    print(f"Found root for x: {root_x}")
                    print()
                break
            except ArithmeticError as err:
                if verbose:
                    print("We could not find a root, " +
                          f"trying other resultants: {err}")
        if root_x is not None:
            break
    if root_x is None:
        print("Cannot find root! Try again with another parameter.")
        return None

    # Solve for the rest,
    eqn_with_x = remove_x_factor(eqn_1.eval_poly(Poly(root_x), x))
    if verbose:
        print("Plugging in x into one of the equations that gave us x:")
        print(eqn_with_x.eval_poly(y, x))
    try:
        root_y = first_root(eqn_with_x, x_guess=root_x)
    except ArithmeticError as err:
        eqn_with_x = remove_x_factor(eqn_2.eval_poly(Poly(root_x), x))
        if verbose:
            print("That failed to converge so we are plugging into the " +
                  "second equation.")
            print(f"Reason: {err}")
            print(eqn_with_x.eval_poly(y, x))
        root_y = first_root(eqn_with_x, x_guess=root_x)
    if verbose:
        print(f"Found root y: {Poly._format_num(root_y)}")
        print("Solving for actual primes next.")
    # Reminder:
    # y = p+q
    # Construct z^2 - 2(p+q)z + pq and our roots are p and q!
    p, q = quadratic(x**2 - x * Poly(root_y) + Poly(n))
    if verbose:
        print(f"Primes recovered: p={
              Poly._format_num(p)} q={Poly._format_num(q)}")
    # Now we can recover d.
    phi = (p-1)*(q-1)
    d = pow(e, -1, phi)
    return p, q, d


def _do_lll(polynomials: list[Poly], x_bound: int, y_bound: int) -> list[Poly]:
    # Create the integer matrix.
    all_monomials = set()
    for poly in polynomials:
        for coeff in poly.all_coeffs():
            all_monomials.add((coeff[1], coeff[2]))
    sorted_monomials = list(all_monomials)
    sorted_monomials.sort(
        key=lambda x: (-1 if x[1] <= x[0] else 1, x[0], x[1]))

    coeffs: list[tuple[int, list[int]]] = []
    for poly in polynomials:
        poly_coeff = [0] * len(sorted_monomials)
        for coeff in poly.all_coeffs():
            col_idx = sorted_monomials.index((coeff[1], coeff[2]))
            scaling_factor = (x_bound ** coeff[1]) * (y_bound ** coeff[2])
            poly_coeff[col_idx] = coeff[0] * scaling_factor
        last_nonzero = len(poly_coeff)-1
        while poly_coeff[last_nonzero] == 0:
            last_nonzero -= 1
        coeffs.append((last_nonzero, poly_coeff))
    coeffs.sort(key=lambda x: x[0])
    L = IntegerMatrix(len(sorted_monomials), len(sorted_monomials))
    for row_idx, (_, poly_coeff) in enumerate(coeffs):
        for col_idx, val in enumerate(poly_coeff):
            L[row_idx, col_idx] = val

    # Run LLL.
    L = LLL.reduction(L)
    # Back to poly.
    result: list[Poly] = []
    for row in L:
        poly = Poly()
        for col_idx, coeff in enumerate(row):
            (x_pow, y_pow) = sorted_monomials[col_idx]
            scaling_factor = (x_bound ** x_pow) * (y_bound ** y_pow)
            if coeff % scaling_factor != 0:
                raise ArithmeticError(
                    "coefficient after LLL is not a multiple of scaling" +
                    "factor")
            poly += Poly(coeff // scaling_factor) * (x ** x_pow) * (y ** y_pow)
        result.append(poly)
    return result


if __name__ == "__main__":
    # Demo: generate a key with small d and attempt Boneh-Durfee recovery
    from key_generation import gen_small_key

    print("[*] Boneh-Durfee Attack Demo")
    print("[*] Generating RSA key with small private exponent...\n")

    e, n, d, phi, p, q = gen_small_key(nbits=128, attempts_per_p=200)
    print(f"n bits: {n.bit_length()}, d bits: {d.bit_length()}")

    print("[*] Running Boneh-Durfee lattice attack (m=4, delta=0.25)...\n")
    res = boneh_durfee_attack(
        n, e, m=3, t=2, delta=0.25, timeout_sec=30, verbose=True)

    if res:
        p_found, q_found, d_found = res
        print("\n[+] SUCCESS! Private key recovered.")
        print(f"    d_actual={d}")
        print(f"    d_found ={d_found}")
        print(f"    Match: {d == d_found}")
    else:
        print("\n[-] Attack did not find the private key.")
        print("    (This may happen for keys outside the vulnerability range.)")
