# ==============================================================================
# Experimental implementation of the Masked SIDH protocol (M-SIDH)
# https://eprint.iacr.org/2023/013.pdf
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

import math
from sage.all import *
from interface import DH_interface, DH_Protocol
from colorama import Back, Style
import logging
import pickle
from sage.misc.persist import SagePickler
import threading
import time
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

proof.all(False)

class MSIDH_Parameters:
    def __init__(self, f, p, E0, A, B, Af, Bf, G, validate=False):
        '''
        Public parameters:
            f: Cofactor for p
            p: A prime such that p = A*B*f - 1
            E0: supersingular curve defined over Fp2
            A, B: two coprime integers that are defines as the sum of t distinct small primes
            Af, Bf: list of factors of A and B (used to speed up calculations down the line)

        Generated parameters:
            PA, QA: Basis of the torsion points of degree A, E0[A] = <PA, QA>
            PB, QB: Basis of the torsion points of degree B, E0[B] = <PB, QB>

        Verification checks:
            A ~ B ~ sqrt(p)
            A, B coprime
            E0 is supersingular
            p is prime and p = A*B*f - 1
            the generated points are on the curve
            the generated points are torsion points of degree A and B
            the generated points are distinct
        '''

        self.f = f
        self.p = p
        self.E0 = E0
        self.A = A
        self.B = B
        self.Af = Af
        self.Bf = Bf
        self.G = G

        # Calculate the points PA, QA, PB, QB

        factorization = factor(p+1)
        print("Factorization of p+1: ", factorization)
        # Sample a random point P on the curve
        P = E0.random_point()

        LP = [ (p+1) / (l**e) * P for l, e in factorization ]
        check_LP = LP.copy()
        for i in range(len(factorization)):
            l, e = factorization[i]
            if e > 1:
                check_LP[i] *= l ** (e-1)
        
        while check_LP.count(E0(0)) != 0:
            print('Restarting from step 1')
            P = E0.random_point()
            LP = [ (p+1) / (l**e) * P for l, e in factorization ]
            check_LP = LP
            for i in range(len(factorization)):
                l, e = factorization[i]
                if e > 1:
                    check_LP[i] *= l ** (e-1)

        # If the order of P is (p+1), then P is a generator of E0
        print(f"Generator P found")

        # SECOND GENERATOR
        Q = E0.random_point()
        while True:
            # Sample a second random point Q on the curve distinct from P
            if Q == P:
                Q = E0.random_point()
                continue

            # Compute the order of Q
            LQ = [ (p+1) / (l**e) * Q for l, e in factorization ]

            check_LQ = LQ.copy()
            for i in range(len(factorization)):
                l, e = factorization[i]
                if e > 1:
                    check_LQ[i] *= l ** (e-1)

            # If the order of Q is not (p+1), then restart from step 7
            if check_LQ.count(E0(0)) != 0:
                Q = E0.random_point()
                print (f"Restarting because LQ.count(E0(0)) != 0, count = {check_LQ.count(E0(0))}")
                continue
            print(f"Found a candidate for Q")
            
            # Check pairwise: multiplicative order of the weil pairing of (LP[i], LQ[i]) is l_i ^ e_i
            #     => Points P,Q are linearly independent
            print("Making P,Q linearly independent...")

            error = False
            for i in range(len(factorization)):
                l, e = factorization[i]
                targets = [l**i for i in range(0, e+1)]
                wp = LP[i].weil_pairing(LQ[i], l**e)
                mo = [wp * i for i in targets]
                if mo.count(1) != 0:
                    print(f"Pair {i} failed, adjusting !")
                    error = True
                    pt_ = (p+1) / (l**e) * E0.random_point()
                    while True:
                        wp_ = pt_.weil_pairing(LP[i], l**e)
                        mo_ = [wp_ * i for i in targets]
                        if  mo_.count(1) != 0:
                            pt_ = (p+1) / (l**e) * E0.random_point()
                            continue
                        else:
                            break
                    Q = Q - pt_

                    
            if not error:
                break
            

        # 13. If the check succeeds, then (P, Q) is a basis of E0[p+1] = <P, Q>
        print(f"Generator Q found")    
        gens = [P, Q]
        self.PA, self.QA = ( B * G * f for G in gens)
        self.PB, self.QB = ( A * G * f for G in gens)

        print(f"{Back.BLUE}==== Generated M-SIDH parameters [{self.__class__.__name__}] ==== {Style.RESET_ALL}")

        # Verify the parameters
        if validate and not self.verify_parameters():
            raise Exception("Invalid parameters")
        

    def __str__(self):
        return f"f: {self.f}\np: {self.p}\nA: {self.A}\nB: {self.B}\nE0: {self.E0}\nPA: {self.PA}\nQA: {self.QA}\nPB: {self.PB}\nQB: {self.QB}"


    def verify_parameters(self):

        p = self.p
        A = self.A
        B = self.B
        f = self.f
        curve = self.E0

        print(f"{Back.LIGHTBLUE_EX}==== Verifying M-SIDH parameters [{self.__class__.__name__}] ==== {Style.RESET_ALL}")
        
        supersingular = curve.is_supersingular(proof=True)
        print(f"Curve is supersingular: {Back.LIGHTGREEN_EX if supersingular else Back.RED}{supersingular}{Style.RESET_ALL}")
        if not supersingular:
            print(f"Curve is not supersingular")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
        # p is prime
        prime = is_prime(p)
        print(f"p is prime: {Back.LIGHTGREEN_EX if prime else Back.RED}{prime}{Style.RESET_ALL}")
        if not prime:
            print(f"p is not prime")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
        # p = A*B*f - 1
    
        valid = p == A*B*f - 1
        print(f"p = A*B*f - 1: {Back.LIGHTGREEN_EX if valid else Back.RED}{valid}{Style.RESET_ALL}")
        if not valid:
            print(f"p != A*B*f - 1")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
        # A ~ B ~ sqrt(p)
        valid = math.sqrt(p)*10e-4 <= A <= math.sqrt(p) * 10e4 and math.sqrt(p) * 10e-4 <= B <= math.sqrt(p) * 10e+4
        print(f"A ~ B ~ sqrt(p): {Back.LIGHTGREEN_EX if valid else Back.RED}{valid}{Style.RESET_ALL}")
        if not valid:
            print(f"A or B is not close to sqrt(p)")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            print(A, B, math.sqrt(p))
            return False
        
        # A, B coprime
        valid = gcd(A, B) == 1
        print(f"A, B coprime: {Back.LIGHTGREEN_EX if valid else Back.RED}{valid}{Style.RESET_ALL}")
        if not valid:
            print(f"A and B are not coprime")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
        PA, QA, PB, QB = self.PA, self.QA, self.PB, self.QB
        
        # points are on the curve
        x, y = self.PA.xy()
        pa_on_curve = curve.is_on_curve(x, y)
        x, y = self.QA.xy()
        qa_on_curve = curve.is_on_curve(x, y)
        x, y = self.PB.xy()
        pb_on_curve = curve.is_on_curve(x, y)
        x, y = self.QB.xy()
        qb_on_curve = curve.is_on_curve(x, y)
        print(f"PA is on the curve: {Back.LIGHTGREEN_EX if pa_on_curve else Back.RED}{pa_on_curve}{Style.RESET_ALL}")
        print(f"QA is on the curve: {Back.LIGHTGREEN_EX if qa_on_curve else Back.RED}{qa_on_curve}{Style.RESET_ALL}")
        print(f"PB is on the curve: {Back.LIGHTGREEN_EX if pb_on_curve else Back.RED}{pb_on_curve}{Style.RESET_ALL}")
        print(f"QB is on the curve: {Back.LIGHTGREEN_EX if qb_on_curve else Back.RED}{qb_on_curve}{Style.RESET_ALL}")
        if not (pa_on_curve and qa_on_curve and pb_on_curve and qb_on_curve):
            print("Some points are not on the curve")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
        # The points have the right order
        valid = PA.order() == A and QA.order() == A and PB.order() == B and QB.order() == B
        print(f"PA has order A: {Back.LIGHTGREEN_EX if valid else Back.RED}{valid}{Style.RESET_ALL}")
        print(f"QA has order A: {Back.LIGHTGREEN_EX if valid else Back.RED}{valid}{Style.RESET_ALL}")
        print(f"PB has order B: {Back.LIGHTGREEN_EX if valid else Back.RED}{valid}{Style.RESET_ALL}")
        print(f"QB has order B: {Back.LIGHTGREEN_EX if valid else Back.RED}{valid}{Style.RESET_ALL}")
        if not valid:
            print("Some points have the wrong order")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            print(PA.order(), QA.order(), PB.order(), QB.order())
            print(A, B)
            return False

        
        # points are distinct
        distinctA = PA != QA
        distinctB = PB != QB
        print(f"PA and QA are distinct: {Back.LIGHTGREEN_EX if distinctA else Back.RED}{distinctA}{Style.RESET_ALL}")
        print(f"PB and QB are distinct: {Back.LIGHTGREEN_EX if distinctB else Back.RED}{distinctB}{Style.RESET_ALL}")
        if not (distinctA and distinctB):
            print("Some points are not distinct")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
        # points are not trivial
        trivial_A = PA == curve(0, 1, 0) or QA == curve(0, 1, 0)
        trivial_B = PB == curve(0, 1, 0) or QB == curve(0, 1, 0)
        print(f"PA and QA are not trivial: {Back.LIGHTGREEN_EX if not trivial_A else Back.RED}{not trivial_A}{Style.RESET_ALL}")
        print(f"PB and QB are not trivial: {Back.LIGHTGREEN_EX if not trivial_B else Back.RED}{not trivial_B}{Style.RESET_ALL}")
        if trivial_A or trivial_B:
            print("Some points are trivial")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
        print(f"{Back.LIGHTGREEN_EX}==== SIDH parameters are valid ==== {Style.RESET_ALL}")
        return True
    

class MSIDHp128(MSIDH_Parameters):
    def __init__(self):
        f =  10
        t = 572
        pari.allocatemem(1<<32)
        print(f"{Back.LIGHTMAGENTA_EX}GENERATING THE SETTINGS...{Style.RESET_ALL}")
        # Get the lambda smallest primes
        primes = Primes()
        primes_list = [primes.unrank(0) ** 2]
        # collect the primes
        for i in range(1,t):
            primes_list.append(primes.unrank(i))

        # A_l = elements of even index in list
        # B_l = elements of odd index in list
        A_l = primes_list[::2]
        B_l = primes_list[1::2]

        # Calculate A and B
        A = prod(A_l)
        B = prod(B_l)

        # Calculate p
        print(f"{Back.LIGHTMAGENTA_EX}CALCULATING THE PRIME...{Style.RESET_ALL}")
        p = A * B * f - 1
        logging.getLogger().setLevel(logging.DEBUG)
        print (f"{Back.LIGHTMAGENTA_EX}GENERATING THE FIELD...{Style.RESET_ALL}")
        F = FiniteField((p, 2), name='x', proof=False)
        print (f"{Back.LIGHTMAGENTA_EX}GENERATING THE CURVE...{Style.RESET_ALL}")
        E0 = EllipticCurve(F, [1,0])
        print(f"{Back.LIGHTMAGENTA_EX}DONE{Style.RESET_ALL}")
        logging.getLogger().setLevel(logging.WARNING)
        super().__init__(f, p, E0, A, B, A_l, B_l, F)


class MSIDHpArbitrary(MSIDH_Parameters):
    def __init__(self, security_parameter, force_t = None):
        self.security_parameter = security_parameter
        t = 2*security_parameter
        if force_t is not None:
            t = force_t
            
        pari.allocatemem(1<<32)
        print(f"{Back.LIGHTMAGENTA_EX}GENERATING THE SETTINGS...{Style.RESET_ALL}")
        # Get the lambda smallest primes
        primes = Primes()
        primes_list = [primes.unrank(0) ** 2]
        # collect the primes
        for i in range(1, 2*t -1):
            primes_list.append(primes.unrank(i))
        # A_l = elements of even index in list
        # B_l = elements of odd index in list
        A_l = primes_list[::2]
        B_l = primes_list[1::2]

        # Calculate A and B
        A = prod(A_l)
        B = prod(B_l)

        # Now we have to find largest n such that prod B <= A_l[n:] ** 2
        n = 0
        while B <= prod(A_l[n:]) ** 2:
            n += 1

        if security_parameter > (t - n + 1):
            # We have to restart with a larger t
            print(f"retrying with t={t+1}")
            self.__init__(security_parameter, force_t=t+1)
            return
        
        # Calculate p
        f = 1
        p = A * B * f - 1
        while not (is_prime(p) and mod(p, 4) == 3):
            f += 1
            p = A * B * f - 1

        assert mod(p, 4) == 3
        print(f"p = {p}")
        F = FiniteField((p, 2), name='x')
        print (f"{Back.LIGHTMAGENTA_EX}GENERATING THE CURVE...{Style.RESET_ALL}")
        (a, ) = F._first_ngens(1)
        print(f"a = {F(1)}")
        E0 = EllipticCurve(j=F(1728))
        print(f"{Back.LIGHTMAGENTA_EX}DONE{Style.RESET_ALL}")
        self.name = f"MSIDH_AES-{security_parameter}"
        super().__init__(f, p, E0, A, B, A_l, B_l, F)


def mewtwo(b, factors):
        '''
        Sample an element x from Z/bZ where x ** 2 = 1 mod b

        '''
        # For each factor, find the square roots
        roots = []
        moduli = []
        for factor in factors:
            # Flip a coin!
            if randrange(2) == 0:
                roots.append(Integer(1))
            else:
                roots.append(Integer(factor - 1))
            moduli.append(Integer(factor))
                
        # Use CRT to find the root
        res = CRT_list(roots, moduli)
        res = IntegerModRing(b)(res)
        assert res ** 2 == 1
        return res
    

class MSIDH_Party_A(DH_interface):
    def __init__(self, parameters):
        self.parameters = parameters

    def get_public_parameters(self):
        return self.parameters
    
    def print_public_parameters(self):
        return f"{self.parameters}"

    def generate_private_key(self):
        alpha = mewtwo(self.parameters.B, self.parameters.Bf)
        a = randrange(self.parameters.A)

        return (alpha, a)

    def compute_public_key(self, private_key):
        pr = self.parameters
        KA = pr.PA + private_key[1] * pr.QA
        phiA = EllipticCurveHom_composite(pr.E0, KA, kernel_order=pr.A)
        return ( phiA.codomain(), private_key[0] * phiA(pr.PB), private_key[0] * phiA(pr.QB) )

    def compute_shared_secret(self, private_key, other_public_key):
        # Check the Weil pairing values
        # eA(Ra, Sa) = eA(PA, QA) ** B
        # Uses sage's implementation of the Weil pairing
        pr = self.parameters
        Ra = other_public_key[1]
        Sa = other_public_key[2]
        p1 = Ra.weil_pairing(Sa, pr.A)
        p2 = pr.PA.weil_pairing(pr.QA, pr.A) ** pr.B
        assert p1 == p2, "Weil pairing values do not match"

        LA = other_public_key[1] + private_key[1] * other_public_key[2]
        psiA = EllipticCurveHom_composite(other_public_key[0], LA, kernel_order=pr.A)
        return psiA.codomain().j_invariant()
    
class MSIDH_Party_B(DH_interface):
    def __init__(self, parameters):
        self.parameters = parameters

    def get_public_parameters(self):
        return self.parameters
    
    def print_public_parameters(self):
        return f"{self.parameters}"

    def generate_private_key(self):
        beta = mewtwo(self.parameters.A, self.parameters.Af)
        b = randrange(self.parameters.B)

        return (beta, b)

    def compute_public_key(self, private_key):
        pr = self.parameters
        KB = pr.PB + private_key[1] * pr.QB
        print("computing isogeny")
        phiB = EllipticCurveHom_composite(pr.E0, KB, kernel_order=pr.B)
        print("Computing public key")
        return ( phiB.codomain(),  private_key[0] * phiB(pr.PA),  private_key[0] * phiB(pr.QA) )

    def compute_shared_secret(self, private_key, other_public_key):
        # Check the Weil pairing values
        # eB(Rb, Sb) = eB(PB, QB) ** A
        # Uses sage's implementation of the Weil pairing
        pr = self.parameters
        Rb = other_public_key[1]
        Sb = other_public_key[2]
        p1 = Rb.weil_pairing(Sb, pr.B)
        p2 = pr.PB.weil_pairing(pr.QB, pr.B) ** pr.A

        assert p1 == p2, "Weil pairing values do not match"

        LB = other_public_key[1] + private_key[1] * other_public_key[2]
        psiB = EllipticCurveHom_composite(other_public_key[0], LB, kernel_order=pr.B)
        return psiB.codomain().j_invariant()


import os.path
def create_protocol(settings_class, additional_parameter=None, p_name=None):
    timer_start = time.time_ns()

    name = settings_class.__name__
    if (not p_name==None) and os.path.isfile(f"{p_name}.pickle"):
        # Load the parameters from file
        print(f"{Back.MAGENTA}Loading {p_name} parameters from file...{Style.RESET_ALL}")
        with open(f"{p_name}.pickle", "rb") as f:
            settings = pickle.load(f)
    else:
        # Generate the parameters
        print(f"{Back.MAGENTA}Generating {name} parameters...{Style.RESET_ALL}")
        if additional_parameter is None:
            settings = settings_class()
        else:
            settings = settings_class(additional_parameter)
        with open(f"{settings.name}.pickle", "wb") as f:
            pickle.dump(settings, f)


    print(f"{Back.GREEN}DONE{Style.RESET_ALL} {(time.time_ns() - timer_start) / 1e9} s")

    partyA = MSIDH_Party_A(settings)
    partyB = MSIDH_Party_B(settings)
    return DH_Protocol(partyA, partyB)


def create_protocol_from_file(path):
    with open(path, "rb") as f:
        settings = pickle.load(f)

    partyA = MSIDH_Party_A(settings)
    partyB = MSIDH_Party_B(settings)
    return DH_Protocol(partyA, partyB)

def create_g128_protocol():
    '''
    Generate the p128 settings
    '''
    time_start = time.time_ns()
    settings = MSIDHp128()
    print(f"{Back.GREEN}DONE{Style.RESET_ALL} {(time.time_ns() - time_start) / 1e9} s")
    with open(f"MSIDHp128.pickle", "wb") as f:
        pickle.dump(settings, f)


def load_128():
    '''
    Load the p128 settings
    '''
    with open(f"MSIDH_AES-8.pickle", "rb") as f:
        settings = pickle.load(f)
    
    return settings