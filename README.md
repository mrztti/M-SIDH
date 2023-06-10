# Masked torsion point Supersingular Isogeny Diffie-Hellman (M-SIDH)
### Implementing countermeasures for attacks on Supersingular Isogeny Diffie-Hellman (SIDH)

**Author:** Malo Ranzetti
**Responsible:** Prof. Serge Vaudenay
**Supervisor:**  Dr. Boris Fouotsa

![LASEC](./figures/logo_lasec_coul.png)

## Abstract
Isogenies between supersingular elliptic curves are useful to construct cryptographic schemes that may be resilient in a post-quantum cryptographic era. One particular scheme proposed is the Supersingular Isogeny Diffie-Hellman (SIDH) key exchange. Since its proposal in 2011, it was seen as a promising candidate. However in 2022 it was shown that one can mount a devastating polynomial time attack against SIDH. Countermeasures to this attack have been proposed by Fuotsa, Moriya and Petit, but imply an explosion in the size of the scheme parameters. 

This project describes the implementation of one countermeasure, namely M-SIDH, for which we implement parameter generation and key exchange for an arbitrary security parameter lambda. After evaluating the system performance, we come to the conclusion that for most practical purposes, these new proposed schemes demand an extreme amount of computational power relative to the security they provide. Making them a viable cryptosystem would require a more efficient algorithm to compute isogenies of large separable degrees. Parameter generation is also affected as we would likely need a third party to pre-compute parameters long in advance.

### [Full project report](./report.pdf)
### [Presentation slides](./presentation.pdf)
### [LASEC homepage](https://lasec.epfl.ch/)

## Usage guide

**Show help:**
    
    sage run.py -h

**Run SIDH implementation on p751 parameters:**
    
    sage run.py -t sidh -c p751

**Generate M-SIDH parameters for lambda = 128 as given in the original M-SIDH paper:**
    
    sage run.py -g128

**Generate M-SIDH parameters for arbitrary lambda given as an argument:**
    
    sage run.py -g <lambda>

**Test 2 rounds of M-SIDH using the parameters for lambda = 128:**
    
    sage run.py -t msidh -r 2 -f MSIDHp128.pickle

**Test 10 rounds of M-SIDH using the parameters for arbitrary lambda = 32:**
    
    sage run.py -t msidh -r 10 -f MSIDH_AES-32.pickle





