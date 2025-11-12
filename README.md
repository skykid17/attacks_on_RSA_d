# SC4010 Project - Attacks on RSA

This project demonstrates Wiener's attack and Boneh Durfee attack on RSA with
small and relatively small private exponents. It includes scripts for generating
both vulnerable and strong RSA key pairs, and for performing the attack to
recover private keys when possible.

## Features
- **Key Generation**: Generate RSA keys that are either vulnerable or resistant
  to Wiener's attack.
- **Wiener's Attack Implementation**: Recover private keys from vulnerable RSA
  moduli using continued fractions.
- **Boneh Durfee Attack Implementation**: Two polynomial time attacks on small
  secret exponent RSA. The attack works when d < N^0.292. The attack is based on
  lattice based Coppersmith's method to find short vectors. It reconstructs
  integer polynomials from the reduced basis (reconstruct_polynomials) and then
  finds small roots using algebraic solvers: pairwise gcds (find_roots_gcd) and
  Groebner‑basis elimination / univariate root finding (find_roots_groebner).
  The quotient trick (pr.quotient(1 + xy - u) and later substituting u = 1 + xy)
  is a linearization used in bivariate constructions (Herrmann & May style) to
  improve speed.
- **Demonstration**: Example runs showing successful and failed attacks.

## Requirements
- PyCryptodome
- fpylll (A)
- SageMath (B)

## Usage A
Install dependencies with:
```
conda create -n venv python=3.12 pycryptodome fpylll
conda activate venv
```
Run the demonstration script:
```
python main.py
```

## Usage B
Install dependencies with:
```
conda config --add channels conda-forge
conda create -n sage_env python=3.12 sage pycryptodome
conda activate sage_env
```
Run the demonstration script:
```
sage -python main.py
```

This will:
- Generate RSA key with small sized private exponent d and attempt Wiener's
  attack and Boneh Durfee attack (both should succeed)
- Generate RSA key with medium sized private exponent d and attempt Wiener's
  attack (should fail) and Boneh Durfee attack (should succeed)
- Generate RSA key with large sized private exponent d and attempt Wiener's
  attack and Boneh Durfee attack (both should fail)

## References
- Wiener, M. J. (1990). Cryptanalysis of Short RSA Secret Exponents. IEEE
  Transactions on Information Theory. https://doi.org/10.1109/18.54902.
- Boneh, D. & Durfee, G. (2000, July). Cryptanalysis of RSA with Private Key d
  Less than N^0.292. IEEE Transactions on Information Theory, 46(4).
  https://doi.org/10.1109/18.850673.
- Herrmann, M. & May, A. (2010, May 26). Maximizing Small Root Bounds by
  Linearization and Applications to Small Secret Exponent RSA.
  https://doi.org/10.1007/978-3-642-13013-7_4.
