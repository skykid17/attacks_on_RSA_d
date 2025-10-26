import time
import math
import secrets
from Crypto.Util.number import getPrime
from typing import Tuple

def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a

def gen_key_with_d_range(
    nbits: int,
    min_d: int,
    max_d: int,
    e_fixed: int = None,
    attempts_per_p: int = 400
) -> Tuple[int,int,int,int,int,int,float]:
    half = nbits // 2
    start = time.perf_counter()

    while True:
        p = getPrime(half)
        for _ in range(attempts_per_p):
            q = getPrime(half)
            if q == p or not (q < p < 2 * q):
                continue
            n = p * q
            phi = (p - 1) * (q - 1)

            try:
                # choose e or respect fixed
                if e_fixed is None:
                    # try small public exponents first
                    for candidate_e in (3, 5, 17, 65537):
                        try:
                            d = pow(candidate_e, -1, phi)
                            e = candidate_e
                            break
                        except ValueError:
                            continue
                    else:
                        continue
                else:
                    e = e_fixed
                    d = pow(e, -1, phi)
            except ValueError:
                continue

            if min_d <= d <= max_d:
                elapsed = time.perf_counter() - start
                return e, n, d, phi, p, q, elapsed

def gen_small_key(nbits: int = 128, attempts_per_p: int = 300, e_fixed: int = None) -> Tuple[int,int,int,int,int,int,float]:
    """
    Generate vulnerable keys (e, n, d, phi, p, q, elapsed_seconds)
    Enforce the classic Wiener bound: d < n^(1/4) / 3
    Returns (e, n, d, phi, p, q, elapsed_seconds)
    """
    half = nbits // 2
    start = time.perf_counter()

    while True:
        p = getPrime(half)
        for _ in range(attempts_per_p):
            q = getPrime(half)
            # ensure p != q and q < p < 2q
            if q == p or not (q < p < 2 * q):
                continue
            n = p * q
            phi = (p - 1) * (q - 1)

            # compute exact maximum allowed d for this n according to classic Wiener's bound:
            # max_d = floor(n^(1/4) / 3)
            max_d = math.isqrt(math.isqrt(n)) // 3
            if max_d < 3:
                continue
            max_d_bits = max(1, max_d.bit_length() - 1)

            # try a few d candidates derived from max_d_bits
            for _d_try in range(12):  # increase a bit the per-(p,q) attempts
                d = secrets.randbits(max_d_bits)
                if d < 3:
                    continue

                # quick filter: ensure d < n^(1/4) / 3 using bit-length heuristic
                if d >= max_d:
                    continue

                if gcd(d, phi) != 1:
                    continue

                # compute e as inverse of d mod phi
                try:
                    e = pow(d, -1, phi)
                except ValueError:
                    continue

                # respect optional fixed e
                if e_fixed is not None and e != e_fixed:
                    continue

                # final exact check for the classic bound
                if d < max_d:
                    elapsed = time.perf_counter() - start
                    return e, n, d, phi, p, q, elapsed

        # loop and generate a new p

def gen_mid_key(nbits: int = 128, attempts_per_p: int = 300) -> Tuple[int,int,int,int,int,int,float]:
    """Generate RSA key with medium d (between Wiener bound and ~n^0.29)."""
    # Compute Wiener's maximum d for a candidate n
    while True:
        p_tmp = getPrime(nbits // 2)
        q_tmp = getPrime(nbits // 2)
        if p_tmp == q_tmp or not (q_tmp < p_tmp < 2 * q_tmp):
            continue
        n_tmp = p_tmp * q_tmp
        wiener_max_d = math.isqrt(math.isqrt(n_tmp)) // 3
        if wiener_max_d < 4:
            continue
        # medium range: (wiener_max_d, n^0.29)
        min_d = wiener_max_d + 1
        max_d = max(min_d + 1, int(n_tmp ** 0.29))
        try:
            return gen_key_with_d_range(nbits=nbits, min_d=min_d, max_d=max_d, attempts_per_p=attempts_per_p)
        except RuntimeError:
            continue

def gen_large_key(nbits: int = 128, e_fixed: int = 65537, attempts_per_p: int = 300) -> Tuple[int,int,int,int,int,int,float]:
    half = nbits // 2
    start = time.perf_counter()

    while True:
        p = getPrime(half)
        q = getPrime(half)
        if p == q or not (q < p < 2 * q):
            continue
        n = p * q
        phi = (p - 1) * (q - 1)

        try:
            d = pow(e_fixed, -1, phi)
        except ValueError:
            continue

        # compute exact minimum allowed d for this n according to classic Wiener's bound:
        # min_d = ceil(n^(1/4) / 3) + 1
        min_d = math.isqrt(math.isqrt(n)) // 3 + 1
        if d > min_d:
            elapsed = time.perf_counter() - start
            return e_fixed, n, d, phi, p, q, elapsed

        # loop and generate a new (p,q) pair

if __name__ == "__main__":
    print("Generating vulnerable key (classic Wiener's bound 81*d^4 < n)...")
    t0 = time.perf_counter()
    e, n, d, phi, p, q, elapsed = gen_small_key(nbits=1024, attempts_per_p=400)
    t_total = time.perf_counter() - t0
    print(f"Done in {t_total:.2f}s (gen elapsed {elapsed:.2f}s).")
    print(f"e: {e}\nd: {d}\nn bits: {n.bit_length()}\nd bits: {d.bit_length()}")
    print("Check Wiener's bound: 81*d^4 < n ->", 81 * pow(d,4) < n)

    print("\nGenerating strong key (d > n^(1/4) / 3)...")
    t0 = time.perf_counter()
    e, n, d, phi, p, q, elapsed = gen_large_key(nbits=1024, e_fixed=65537)
    t_total = time.perf_counter() - t0
    print(f"Done in {t_total:.2f}s (gen elapsed {elapsed:.2f}s).")
    print(f"e: {e}\nd: {d}\nn bits: {n.bit_length()}\nd bits: {d.bit_length()}")
    print("Check Wiener's bound: 81*d^4 < n ->", 81 * pow(d,4) < n)

