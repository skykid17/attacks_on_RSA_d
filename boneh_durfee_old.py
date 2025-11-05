from typing import Optional, Tuple, List
import math
import time
from key_generation import gcd
from fpylll import IntegerMatrix, LLL
from fpylll.gso import MatGSO

def _is_perfect_square(n: int) -> bool:
    """Return True if n is a perfect square."""
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n

def boneh_durfee_attack(
    n: int,
    e: int,
    m: int = 4,
    X_bound: Optional[int] = None,
    Y_bound: Optional[int] = None,
    delta: float = 0.18,
    timeout_sec: float = 60,
) -> Optional[Tuple[int, int, int]]:
    if n % 2 == 0:
        # n must be even, as such, we can already factorize it!
        phi = n // 2 - 1
        return n // 2, 2, pow(e, -1, phi)

    # Some default values. Y should be in the range p+q and X should be in the
    # range e^delta.
    if X_bound is None:
        X_bound = 2 ** (int(n.bit_length() * delta) + 1)
    if Y_bound is None:
        Y_bound = 2 ** (n.bit_length() // 2 + 1)

    # We are trying to solve:
    # f(x, y) = x(A-y)-1 = 0 (mod e).
    # A = N+1, y = p+q
    #
    # This is derived from:
    # ed = -k phi(n) + 1
    # ed = -k(N-p-q+1) + 1
    # k(N-p-q+1) = 1 (mod e).
    A = n + 1
    big_x = x * Poly(X_bound)
    big_y = y * Poly(Y_bound)
    f = big_x * (Poly(A) - big_y) - Poly(1)
    e_poly = Poly(e)

    polynomials: list[Poly] = []
    # Boneh and Durfee found that the best m is 7 and t is 3.
    # That means our equations are 0 (mod e^7) and have a maximum power of y^3.
    t = 3
    for k in range(m + 1):
        # x-shifts: x^i * f(x, y)^k * e^(m-k) = 0 (mod e^m)
        # For each k, we use g_ik for i = 0, ..., m-k.
        for i in range(m-k+1):
            polynomials.append(((big_x**i) * (f**k) * (e_poly ** (m-k))))
        # y-shifts: y^i * f(x, y)^k * e^(m-k) = 0 (mod e^m)
        # For each k, we use h_lk for l = 1, ..., t.
        for l_pow in range(1, t+1):
            polynomials.append(((big_y**l_pow) * (f**k) * (e_poly ** (m-k))))

    polynomials = _do_lll(polynomials)
    print("Equations after LLL:")
    eqn_1, eqn_2, eqn_3 = polynomials[1], polynomials[2], polynomials[3]
    print(eqn_1)
    print(eqn_2)
    print(eqn_3)
    print("Resultant was: ")
    print(remove_x_factor(eqn_1.resultant(eqn_2)))
    print(remove_x_factor(eqn_1.resultant(eqn_3)))
    print(remove_x_factor(eqn_2.resultant(eqn_3)))
    return None


def _do_lll(polynomials: list[Poly]) -> list[Poly]:
    # Create the integer matrix.
    all_monomials = set()
    for poly in polynomials:
        for coeff in poly.all_coeffs():
            all_monomials.add((coeff[1], coeff[2]))
    sorted_monomials = list(all_monomials)
    sorted_monomials.sort(
        key=lambda x: (x[0]+x[1], x[0], x[1]))

    L = IntegerMatrix(len(sorted_monomials), len(sorted_monomials))
    for row_idx, poly in enumerate(polynomials):
        for coeff in poly.all_coeffs():
            col_idx = sorted_monomials.index((coeff[1], coeff[2]))
            L[row_idx, col_idx] = coeff[0]

    # Run LLL.
    L = LLL.reduction(L)
    # Back to poly.
    result: list[Poly] = []
    for row in L:
        poly = Poly()
        for col_idx, coeff in enumerate(row):
            (x_pow, y_pow) = sorted_monomials[col_idx]
            poly += Poly(coeff) * (x ** x_pow) * (y ** y_pow)
        result.append(poly)
    return result


if __name__ == "__main__":
    # Demo: generate a key with small d and attempt Boneh-Durfee recovery
    from key_generation import gen_small_key

    print("[*] Boneh-Durfee Attack Demo")
    print("[*] Generating RSA key with small private exponent...\n")

    e, n, d, phi, p, q = gen_small_key(nbits=16, attempts_per_p=200)
    print(f"n bits: {n.bit_length()}, d bits: {d.bit_length()}")

    print("[*] Running Boneh-Durfee lattice attack (m=4, delta=0.18)...\n")
    res = boneh_durfee_attack(n, e, m=4, delta=0.18, timeout_sec=30)

    if res:
        p_found, q_found, d_found = res
        print("\n[+] SUCCESS! Private key recovered.")
        print(f"    d_actual={d}")
        print(f"    d_found ={d_found}")
        print(f"    Match: {d == d_found}")
    else:
        print("\n[-] Attack did not find the private key.")
        print("    (This may happen for keys outside the vulnerability range.)")
