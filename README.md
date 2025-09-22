# Wiener's Attack on RSA - SC4010 Project

This project demonstrates Wiener's attack on RSA cryptosystems with small private exponents. It includes scripts for generating both vulnerable and strong RSA key pairs, and for performing the attack to recover private keys when possible.

## Features
- **Key Generation**: Generate RSA keys that are either vulnerable or resistant to Wiener's attack.
- **Wiener's Attack Implementation**: Recover private keys from vulnerable RSA moduli using continued fractions.
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
