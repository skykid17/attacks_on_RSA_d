import time
import math
import secrets
from Crypto.Util.number import getPrime
from typing import Tuple

def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a

def gen_vulnerable_keys_fast(
    nbits: int = 1024,
    attempts_per_p: int = 400,
    e_fixed: int = None
) -> Tuple[int,int,int,int,int,int,float]:
    """
    Generate vulnerable keys (e, n, d, phi, p, q, elapsed_seconds)
    Enforce the classic Wiener bound: d < n^(1/4) / 3
    - nbits: target bitlength of modulus n
    - attempts_per_p: number of q tries per p
    - e_fixed: optional fixed public exponent
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

if __name__ == "__main__":
    print("Generating vulnerable key (classic Wiener's bound 81*d^4 < n)...")
    t0 = time.perf_counter()
    e, n, d, phi, p, q, elapsed = gen_vulnerable_keys_fast(nbits=1024, attempts_per_p=400)
    t_total = time.perf_counter() - t0
    print(f"Done in {t_total:.2f}s (gen elapsed {elapsed:.2f}s).")
    print(f"e: {e}\nd: {d}\nn bits: {n.bit_length()}\nd bits: {d.bit_length()}")
    print("Check Wiener's bound: 81*d^4 < n ->", 81 * pow(d,4) < n)
