import math
import secrets
import time
from decimal import Decimal
from typing import Tuple

from Crypto.Util.number import getPrime, getRandomRange


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def gen_p_and_q(nbits: int, attempts_per_p: int = 400) -> Tuple[int, int]:
    half = nbits // 2

    while True:
        p = getPrime(half)
        if p.bit_length() < half:
            continue
        for _ in range(attempts_per_p):
            q = getPrime(half)
            if q == p or q.bit_length() < half:
                continue
            n = p * q
            if n.bit_length() < nbits:
                continue
            return p, q


def gen_key_with_d_range(
    nbits: int, min_d: int, max_d: int, e_fixed: int = None, attempts_per_p: int = 400
) -> Tuple[int, int, int, int, int, int]:
    while True:
        p, q = gen_p_and_q(nbits, attempts_per_p)
        n = p * q
        if n.bit_length() < nbits:
            continue
        phi = (p - 1) * (q - 1)

        try:
            candidate_d = getRandomRange(min_d, max_d)
            e = pow(candidate_d, -1, phi)
            return e, n, candidate_d, phi, p, q
        except ValueError:
            continue


def gen_small_key(
    nbits: int = 128, attempts_per_p: int = 300, e_fixed: int = None
) -> Tuple[int, int, int, int, int, int]:
    """
    Generate vulnerable keys (e, n, d, phi, p, q)
    Enforce the classic Wiener bound: d < n^(1/4) / 3
    Returns (e, n, d, phi, p, q)
    """
    max_d = math.isqrt(math.isqrt(2**nbits)) // 3
    return gen_key_with_d_range(nbits=nbits, min_d=2, max_d=max_d, attempts_per_p=attempts_per_p)


def gen_mid_key(nbits: int = 128, attempts_per_p: int = 300) -> Tuple[int, int, int, int, int, int]:
    """Generate RSA key with medium d (between Wiener bound and ~n^0.29)."""
    # Compute Wiener's maximum d for a candidate n
    min_d = math.isqrt(math.isqrt(2**nbits)) // 3
    max_d = int(Decimal(2 ** (nbits - 1)) ** Decimal(0.28))
    return gen_key_with_d_range(
        nbits=nbits, min_d=min_d, max_d=max_d, attempts_per_p=attempts_per_p
    )


def gen_large_key(
    nbits: int = 128, e_fixed: int = 65537, attempts_per_p: int = 300
) -> Tuple[int, int, int, int, int, int]:
    while True:
        p, q = gen_p_and_q(nbits, attempts_per_p)
        n = p * q
        phi = (p - 1) * (q - 1)
        try:
            d = pow(e_fixed, -1, phi)
        except ValueError:
            continue

        return e_fixed, n, d, phi, p, q


if __name__ == "__main__":
    print("Generating vulnerable key (classic Wiener's bound 81*d^4 < n)...")
    t0 = time.perf_counter()
    e, n, d, phi, p, q, elapsed = gen_small_key(nbits=1024, attempts_per_p=400)
    t_total = time.perf_counter() - t0
    print(f"Done in {t_total:.2f}s (gen elapsed {elapsed:.2f}s).")
    print(f"e: {e}\nd: {d}\nn bits: {
          n.bit_length()}\nd bits: {d.bit_length()}")
    print("Check Wiener's bound: 81*d^4 < n ->", 81 * pow(d, 4) < n)

    print("\nGenerating strong key (d > n^(1/4) / 3)...")
    t0 = time.perf_counter()
    e, n, d, phi, p, q = gen_large_key(nbits=1024, e_fixed=65537)
    t_total = time.perf_counter() - t0
    print(f"Done in {t_total:.2f}s (gen elapsed {elapsed:.2f}s).")
    print(f"e: {e}\nd: {d}\nn bits: {
          n.bit_length()}\nd bits: {d.bit_length()}")
    print("Check Wiener's bound: 81*d^4 < n ->", 81 * pow(d, 4) < n)
