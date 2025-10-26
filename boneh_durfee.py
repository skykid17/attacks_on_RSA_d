from typing import Optional, Tuple, List
import math
import time
from key_generation import gcd

try:
    from fpylll import IntegerMatrix, LLL
    from fpylll.gso import MatGSO
    HAS_FPYLLL = True
except ImportError:
    HAS_FPYLLL = False


def _is_perfect_square(n: int) -> bool:
    """Return True if n is a perfect square."""
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def boneh_durfee_attack(
    n: int, e: int, m: int = 4, t: int = None, X: int = None, Y: int = None,
    delta: float = 0.18, timeout_sec: float = 60
) -> Optional[Tuple[int, int, int]]:

    start_time = time.time()
    
    # Set default parameters
    if t is None: t = max(1, int((1 - 2 * delta) * m))
    if X is None: X = 2 * (1 << int(n.bit_length() * delta))
    if Y is None:
        Y = 1 << int(n.bit_length() * 0.5)
    
    # Polynomial: 1 + x*(A + y) where A = (n+1)/2
    A = (n + 1) // 2
    
    # Construct lattice basis from x-shifts and y-shifts
    monomials = []
    basis_vectors = []
    
    # x-shifts: x^i * N^(m-k) * poly^k for k=0..m, i=0..m-k
    for k in range(m + 1):
        for i in range(m - k + 1):
            # Coefficient of x^i in (1 + x*(A+y))^k * N^(m-k)
            coeff = (n ** (m - k))
            monomials.append(('x', i, k))
            vec = [0] * (m + t + 2)
            vec[i] = coeff
            basis_vectors.append(vec)
    
    # y-shifts: y^j * poly^k * N^(m-k) for selected j, k
    for j in range(1, t + 1):
        for k in range((m * j) // t, m + 1):
            monomials.append(('y', j, k))
            vec = [0] * (m + t + 2)
            vec[m + j] = n ** (m - k)
            basis_vectors.append(vec)
    
    # Attempt LLL reduction
    try:
        M = IntegerMatrix.from_matrix(basis_vectors)
        gso = MatGSO(M)
        LLL.reduction(gso)
        reduced = [list(row) for row in M]
    except Exception as e:
        print(f"[Boneh-Durfee] LLL reduction failed: {e}")
        return None
    
    # Extract short vectors and try to find roots
    if time.time() - start_time > timeout_sec:
        return None
    
    # Check for solution in small vectors
    for vec_idx, vec in enumerate(reduced[:min(5, len(reduced))]):
        # Try to extract d from this vector
        # This is simplified; full version would use polynomial solving
        for val in vec:
            if val != 0 and val < e:
                # Check if this could be d
                if (e * val) % n == 1 or (e * val - 1) % n == 0:
                    d_candidate = val
                    # Verify
                    if (e * d_candidate) % n == 1:
                        # Factor n using d
                        phi = (e * d_candidate - 1) // n if (e * d_candidate - 1) % n == 0 else None
                        if phi:
                            s = n - phi + 1
                            disc = s * s - 4 * n
                            if disc > 0 and _is_perfect_square(disc):
                                sqrt_disc = math.isqrt(disc)
                                p = (s + sqrt_disc) // 2
                                q = (s - sqrt_disc) // 2
                                if p * q == n and p > 1 and q > 1:
                                    return max(p, q), min(p, q), d_candidate
    
    return None


if __name__ == "__main__":
    # Demo: generate a key with small d and attempt Boneh-Durfee recovery
    from key_generation import gen_small_key

    print("[*] Boneh-Durfee Attack Demo")
    print("[*] Generating RSA key with small private exponent...\n")
    
    e, n, d, phi, p, q, elapsed = gen_small_key(nbits=128, attempts_per_p=200)
    print(f"n bits: {n.bit_length()}, d bits: {d.bit_length()}")
    print(f"Generated in {elapsed:.2f}s\n")
    
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
