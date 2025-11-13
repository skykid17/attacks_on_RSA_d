#!/usr/bin/python3
from math import isqrt


def cf_expansion(n, d):
    e = []
    while d:
        q = n // d
        e.append(q)
        n, d = d, n % d
    return e


def convergents(cf):
    h1, h2 = 1, 0
    k1, k2 = 0, 1
    for a in cf:
        h = a * h1 + h2
        k = a * k1 + k2
        yield h, k
        h2, h1 = h1, h
        k2, k1 = k1, k


def wiener_attack(n, e):
    cf = cf_expansion(e, n)
    for k, d in convergents(cf):
        if k == 0:
            continue
        if (e * d - 1) % k == 0:
            phi = (e * d - 1) // k
            s = n - phi + 1
            discriminant = s * s - 4 * n
            if discriminant >= 0:
                sqrt_disc = isqrt(discriminant)
                if sqrt_disc * sqrt_disc == discriminant:
                    p = (s + sqrt_disc) // 2
                    q = (s - sqrt_disc) // 2
                    if p * q == n:
                        return max(p, q), min(p, q), d
    return None


if __name__ == "__main__":
    # Lightweight self-check to ensure the continued fraction logic holds up.
    import time
    import sys
    from pathlib import Path
    if __package__ in {None, ""}:
        sys.path.append(str(Path(__file__).resolve().parents[1]))
    from mathlib.key_generation import gen_large_key, gen_mid_key, gen_small_key
    
    demo_start = time.perf_counter()
    e, n, d, _phi, p, q = gen_small_key(nbits=2048, attempts_per_p=80)
    print(
        f"Generated RSA modulus with n_bits={n.bit_length()} d_bits={d.bit_length()}"
    )
    print("Running Wiener attack...\n")

    attack_start = time.perf_counter()
    res = wiener_attack(n, e)
    attack_time = time.perf_counter() - attack_start

    if res is None:
        print(f"Attack failed after {attack_time:.2f}s. Consider regenerating the key.")
    else:
        p_found, q_found, d_found = res
        print("\n[+] SUCCESS! Private key recovered.")
        print(f"    d_actual={d}")
        print(f"    d_found ={d_found}")
        print(f"    Match: {d == d_found}")

    total_time = time.perf_counter() - demo_start
    print(f"Demo finished in {total_time:.2f}s.")