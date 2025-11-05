#!/usr/bin/env python
"""Main runner comparing Wiener's attack with Boneh-Durfee attack.

This demonstrates RSA key vulnerability analysis:
- Small d (Wiener-vulnerable): Wiener's attack should succeed
- Medium d: Wiener fails, Boneh-Durfee may find with bounded search
- Large d (strong): Both attacks fail

Boneh-Durfee uses pure-Python polynomial root-finding (no external LLL libraries).
"""

import time
from typing import Optional, Tuple

from boneh_durfee import boneh_durfee_attack
from key_generation import gen_large_key, gen_mid_key, gen_small_key
from wiener_attack import wiener_attack


def run_case(
    e: int, n: int, d: int
) -> Tuple[bool, bool, Optional[Tuple[int, int, int]], Optional[Tuple[int, int, int]]]:
    """Run Wiener and Boneh-Durfee on given public key and return results.

    Returns a tuple: (wiener_success, bd_success, wiener_result, bd_result)
    where results are (p,q,d) or None.
    """
    print(f"n bits: {n.bit_length()}, d bits: {d.bit_length()}\n")

    # Wiener
    t0 = time.perf_counter()
    wres = wiener_attack(n, e)
    t_w = time.perf_counter() - t0
    w_ok = wres is not None
    print(f"Wiener: {'succeeded' if w_ok else 'failed'} (time {t_w:.3f}s)")
    if w_ok:
        wp, wq, wd = wres
        print(f"d={wd}\np bits={wp.bit_length()} q bits={wq.bit_length()}\n")

    # Boneh-Durfee
    # For keys with small d, k is large, so we need larger search space.
    # For demonstration: use moderate max_k
    t0 = time.perf_counter()
    # exit after 86400 seconds (24 hours)
    try:
        bd_res = boneh_durfee_attack(n, e, n.bit_length()//2, m=3, delta=0.26)
    except TimeoutError:
        print("Boneh-Durfee: failed (timeout)")
        bd_ok = False
        return w_ok, bd_ok, wres, None
    t_bd = time.perf_counter() - t0
    bd_ok = bd_res is not None
    print(f"Boneh-Durfee: {'succeeded' if bd_ok else 'failed'} (time {t_bd:.3f}s)")
    if bd_ok:
        bp, bq, bd = bd_res
        print(f"d={bd}\np bits={bp.bit_length()} q bits={bq.bit_length()}\n")

    return w_ok, bd_ok, wres, bd_res


def main():
    """Run comparative attack demonstrations on three key types."""
    nbits = 1024

    # 1) Small d (Wiener vulnerable)
    print("\n=== Case [1/3]: Small d ===")
    e_s, n_s, d_s, phi_s, p_s, q_s = gen_small_key(nbits=nbits)
    w_s, bd_s, wres_s, bdres_s = run_case(e_s, n_s, d_s)

    # 2) Medium d (Wiener should fail, Boneh-Durfee succeed)
    print("\n=== Case [2/3]: Medium d ===")
    e_m, n_m, d_m, phi_m, p_m, q_m = gen_mid_key(nbits=nbits)
    w_m, bd_m, wres_m, bdres_m = run_case(e_m, n_m, d_m)

    # 3) Large d (strong key, both fail)
    print("\n=== Case [3/3]: Large d ===")
    e_l, n_l, d_l, phi_l, p_l, q_l = gen_large_key(nbits=nbits)
    w_l, bd_l, wres_l, bdres_l = run_case(e_l, n_l, d_l)


if __name__ == "__main__":
    main()
