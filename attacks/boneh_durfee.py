from decimal import Decimal
from typing import Optional, Tuple
from fpylll import LLL, IntegerMatrix
from mathlib.utils import Poly, x, y, z, remove_x_factor, resultant, first_root, quadratic

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
    X_assert: Optional[int] = None,
    Y_assert: Optional[int] = None,
) -> Optional[Tuple[int, int, int]]:
    if verbose:
        print(f"Solving n={Poly._format_num(n)} e={Poly._format_num(e)}")
    assert n % 2 != 0
    # Some default values. Y should be in the range p+q and X should be in the
    # range e^delta.
    if X_bound is None:
        X_bound = int(Decimal(n) ** Decimal(delta))
        if verbose:
            print("X_bound:", X_bound.bit_length())
    if Y_bound is None:
        # Y_bound = 4 * int(Decimal(n) ** Decimal(0.5))
        Y_bound = 4 * int(Decimal(n) ** Decimal(0.5))
        if verbose:
            print("Y_bound:", Y_bound.bit_length())
    U_bound = X_bound * Y_bound + 1
    if verbose:
        print("U_bound:", U_bound.bit_length())

    # We are trying to solve:
    # f(x, y) = x(A+y)+1 = 0 (mod e).
    #         = Ax + u
    # A = N+1, y = -p-q, u = xy+1
    #
    # This is derived from:
    # ed = k phi(n) + 1
    # ed = k(N-p-q+1) + 1
    # k(N-p-q+1) + 1 = 0 (mod e).
    A = n + 1
    u = x * y + 1
    f = A * x + z
    if X_assert is not None and Y_assert is not None:
        U_assert = X_assert*Y_assert+1
        assert f.eval(X_assert, Y_assert, U_assert) % e == 0

    if verbose:
        # Use a simple formatted string to avoid multi-line f-string parsing
        # differences between Python interpreters.
        print("Function we are trying to solve is (mod e={0}):".format(
            Poly._format_num(e)))
        print(f)
        print()

    polynomials: list[Poly] = []
    # Boneh and Durfee found that the best m is 7 and t is 3.
    # That means our equations are 0 (mod e^7) and have a maximum power of y^3.
    f_pow: list[int | Poly] = [1, f]
    fp = f
    for _ in range(2, m+1):
        fp *= f
        f_pow.append(fp)
    for k in range(m + 1):
        # x-shifts: x^i * f(x, y, z)^k * e^(m-k) = 0 (mod e^m)
        # For each k, we use g_ik for i = 0, ..., m-k.
        for i in range(m-k+1):
            g_ik = (x**i) * f_pow[k] * (e ** (m-k))
            polynomials.append(g_ik)
            if X_assert is not None and Y_assert is not None:
                assert g_ik.eval(X_assert, Y_assert, U_assert) % (e ** m) == 0
    for k in range(m + 1):
        # y-shifts: y^i * f(x, y, z)^k * e^(m-k) = 0 (mod e^m)
        # For each k, we use h_lk for l = 1, ..., t.
        for l_pow in range(1, t+1):
            h_lk = (y**l_pow) * f_pow[k] * (e ** (m-k))
            polynomials.append(h_lk)
            if X_assert is not None and Y_assert is not None:
                assert h_lk.eval(X_assert, Y_assert, U_assert) % (e ** m) == 0
    polynomials = [poly.sub_xy(z - 1) for poly in polynomials]
    if X_assert is not None and Y_assert is not None:
        for poly in polynomials:
            assert poly.eval(X_assert, Y_assert, U_assert) % (e ** m) == 0

    if verbose:
        print("Equations before LLL:")
        for i in range(4):
            print(polynomials[i])
        print()
        print("Doing LLL...")
        print()
    polynomials = _do_lll(polynomials, X_bound, Y_bound, U_bound, verbose)
    if verbose:
        print("Equations after LLL:")
        for i in range(4):
            print(polynomials[i])
        print()
    if X_assert is not None and Y_assert is not None:
        num_valid = sum([
            1 if poly.eval(X_assert, Y_assert, U_assert) == 0
            else 0
            for poly in polynomials
        ])
        print(f"Number of valid polys: {num_valid}, they are ")
        for i, poly in enumerate(polynomials):
            if poly.eval(X_assert, Y_assert, U_assert) == 0:
                print(i, end=" ")
        print()
        assert num_valid >= 2
    root_x = None
    for i in range(1, 5):
        if verbose:
            print(f"Trying {i=}")
        poly_i = polynomials[i].eval_poly(x, y, u)
        for j in range(i):
            poly_j = polynomials[j].eval_poly(x, y, u)
            r_candidate = resultant(poly_i, poly_j)
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
                if X_assert is not None and Y_assert is not None:
                    if (eqn_1.eval(X_assert, Y_assert, 0) == 0 and
                            eqn_2.eval(X_assert, Y_assert, 0) == 0):
                        assert r_candidate.eval(X_assert, Y_assert, 0) == 0
                root_x = first_root(
                    remove_x_factor(r_candidate),
                    x_guess=X_bound//2)
                if X_assert is not None:
                    assert X_assert == root_x
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
    eqn_with_x = remove_x_factor(eqn_1.eval_poly(root_x, -x, x))
    if verbose:
        print("Plugging in x into one of the equations that gave us x:")
        print(eqn_with_x.eval_poly(y, x, x))
    if Y_assert is not None:
        assert eqn_with_x.eval(-Y_assert, 0, 0) == 0
    try:
        root_y = first_root(eqn_with_x, x_guess=Y_bound)
    except ArithmeticError as err:
        eqn_with_x = remove_x_factor(eqn_2.eval_poly(root_x, -x, x))
        if verbose:
            print("That failed to converge so we are plugging into the " +
                  "second equation.")
            print(f"Reason: {err}")
            print(eqn_with_x.eval_poly(y, x, x))
        if Y_assert is not None:
            assert eqn_with_x.eval(-Y_assert, 0, 0) == 0
        root_y = first_root(eqn_with_x, x_guess=Y_bound)
    if verbose:
        print(f"Found root y: {Poly._format_num(root_y)}")
        print("Solving for actual primes next.")
    # Reminder:
    # y = -p-q
    # Construct z^2 - (p+q)z + pq and our roots are p and q!
    p, q = quadratic(x**2 - x * root_y + n)
    if verbose:
        print("Primes recovered: p={0} q={1}".format(
            Poly._format_num(p), Poly._format_num(q)))
    # Now we can recover d.
    phi = (p-1)*(q-1)
    d = pow(e, -1, phi)
    return p, q, d


def _do_lll(
        polynomials: list[Poly],
        x_bound: int, y_bound: int, u_bound: int, verbose: bool) -> list[Poly]:
    # Create the integer matrix.
    all_monomials = set()
    for poly in polynomials:
        for (_, pow) in poly.all_coeffs():
            all_monomials.add(pow)
    sorted_monomials = list(all_monomials)
    sorted_monomials.sort(key=lambda x: (
        # x-shifts first. Anything with y goes to the back.
        1 if x[1] > 0 else -1,
        x[0]+x[2],  # Sum of x power and u power.
        x[2], x[0],  # Smaller u power first (x before u. x^2 before ux before
                     # u^2).
        x[1],  # Finally, y power.
    ))

    coeffs: list[tuple[int, list[int]]] = []
    for poly in polynomials:
        poly_coeff = [0] * len(sorted_monomials)
        for (coeff, pow) in poly.all_coeffs():
            col_idx = sorted_monomials.index(pow)
            scaling_factor = (
                (x_bound ** pow[0]) *
                (y_bound ** pow[1]) *
                (u_bound ** pow[2])
            )
            poly_coeff[col_idx] = coeff * scaling_factor
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
            pow = sorted_monomials[col_idx]
            (x_pow, y_pow, u_pow) = pow
            scaling_factor = (
                (x_bound ** x_pow) *
                (y_bound ** y_pow) *
                (u_bound ** u_pow)
            )
            if coeff % scaling_factor != 0:
                raise ArithmeticError(
                    "coefficient after LLL is not a multiple of scaling" +
                    "factor")
            poly[pow] += coeff // scaling_factor
        result.append(poly)
    return result


if __name__ == "__main__":
    # Demo: generate a key with small d and attempt Boneh-Durfee recovery
    import time
    from mathlib.key_generation import gen_large_key, gen_mid_key, gen_small_key

    demo_start = time.perf_counter()
    print("[*] Boneh-Durfee Attack Demo")
    print("[*] Generating RSA key with small private exponent...\n")

    e, n, d, phi, p, q = gen_mid_key(nbits=1024, attempts_per_p=200)
    print(f"n bits: {n.bit_length()}, d bits: {d.bit_length()}")

    actual_x = (e*d - 1) // (p-1) // (q-1)
    actual_y = -p-q

    print("[*] Running Boneh-Durfee lattice attack (m=4, delta=0.26)...\n")
    res = boneh_durfee_attack(
        n, e, m=4, t=3, delta=0.25, timeout_sec=30, verbose=True,
        X_assert=actual_x, Y_assert=actual_y)

    if res:
        p_found, q_found, d_found = res
        print("\n[+] SUCCESS! Private key recovered.")
        print(f"    d_actual={d}")
        print(f"    d_found ={d_found}")
        print(f"    Match: {d == d_found}")
    else:
        print("\n[-] Attack did not find the private key.")
        print("    (This may happen for keys outside the vulnerability range.)")
    total_time = time.perf_counter() - demo_start
    print(f"\nDemo finished in {total_time:.2f}s.")