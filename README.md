# SC4010 Project - Attacks on RSA

This project demonstrates Wiener's attack, as well as Boneh Durfee on RSA cryptosystems with small and relatively small private exponents. It includes scripts for generating both vulnerable and strong RSA key pairs, and for performing the attack to recover private keys when possible.

## Features
- **Key Generation**: Generate RSA keys that are either vulnerable or resistant to Wiener's attack.
- **Wiener's Attack Implementation**: Recover private keys from vulnerable RSA moduli using continued fractions.
- **Boneh Durfee Attack Implementation**: Two polynomial time attacks on small secret exponent RSA. The attack works when d < N^0.292. The attack is based on lattice based Coppersmith's method to solve modular equations.

- **Demonstration**: Example runs showing successful and failed attacks.

## Requirements
Install dependencies with:
```
pip install -r requirements.txt
```

## Usage
Run the demonstration script:
```
python wiener_attack.py
```
This will:
- Generate a vulnerable RSA key and attempt Wiener's attack (should succeed)
- Generate a strong RSA key and attempt Wiener's attack (should fail)

## References
- M. Wiener, "Cryptanalysis of Short RSA Secret Exponents," IEEE Transactions on Information Theory, 1990.

---
