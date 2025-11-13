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
- **Boneh-Durfee Attack Implementation**: Single streamlined pipeline that
  builds the Herrmann–May style lattice, reduces it with LLL, and then uses
  resultant elimination plus an Aberth–Ehrlich root search (all in pure Python)
  to recover the short secret exponent when `d < N^0.29`.
- **Demonstration**: Example runs showing successful and failed attacks.

## Requirements
- PyCryptodome
- fpylll

## Usage
Install dependencies with:
```
conda create -n venv python=3.12 pycryptodome fpylll
conda activate venv
```
Run the demonstration script:
```
python main.py
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
- Aberth, O. (1973). Iteration methods for finding all zeros of a polynomial
  simultaneously. Matematički Vesnik, 10(1). (Aberth–Ehrlich method used for
  integer root polishing in `mathlib`).
