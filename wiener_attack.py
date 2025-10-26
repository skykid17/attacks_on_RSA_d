#!/usr/bin/python3
from math import isqrt
from key_generation import gen_small_key, gen_large_key

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

if __name__ == '__main__':
    # Generate vulnerable RSA key
    e, n, d, phi, p, q, elapsed = gen_small_key(nbits=2048, attempts_per_p=400)
    print(f"Generated vulnerable RSA key in {elapsed:.2f} seconds:")
    print(f"p = {p}\nq = {q}\nd = {d}\n")
    result = wiener_attack(n, e)
    if result:
        print("Wiener's attack successful!")
        p, q, d = result
        print(f"Found factors: \np = {p}\nq = {q}\nd = {d}")
    else:
        print("Wiener's attack failed; could not factor n.")
    e, n, d, phi, p, q, elapsed = gen_large_key(nbits=2048, e_fixed=65537)
    print(f"Generated strong RSA key in {elapsed:.2f} seconds:")
    print(f"p = {p}\nq = {q}\nd = {d}\n")
    result = wiener_attack(n, e)
    if result:
        print("Wiener's attack unexpectedly succeeded on strong key!")
        p, q, d = result
        print(f"Found factors: \np = {p}\nq = {q}\nd = {d}")
    else:
        print("Wiener's attack failed on strong key; could not factor n.")