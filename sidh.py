# ==============================================================================
# Implementation of the Supersingular Isogeny Diffie-Hellman key exchange
# WARNING: Vulnerable to devastating polynomial time attack
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

from sage.all import *
from interface import DH_interface, DH_Protocol
from colorama import Back, Style

class SIDH_Party_A(DH_interface):
    def __init__(self, parameters):
        self.parameters = parameters

    def get_public_parameters(self):
        return self.parameters
    
    def print_public_parameters(self):
        return f"{self.parameters}"

    def generate_private_key(self):
        return randrange(self.parameters.lA ^ self.parameters.eA)

    def compute_public_key(self, private_key):
        pr = self.parameters
        KA = pr.PA + private_key * pr.QA
        phiA = pr.curve.isogeny(KA, algorithm="factored")
        return ( phiA.codomain(), phiA(pr.PB), phiA(pr.QB) )

    def compute_shared_secret(self, private_key, other_public_key):
        LA = other_public_key[1] + private_key * other_public_key[2]
        psiA = other_public_key[0].isogeny(LA, algorithm="factored")
        return psiA.codomain().j_invariant()
    
class SIDH_Party_B(DH_interface):
    def __init__(self, parameters):
        self.parameters = parameters

    def get_public_parameters(self):
        return self.parameters
    
    def print_public_parameters(self):
        return f"{self.parameters}"

    def generate_private_key(self):
        return randrange(self.parameters.lB ^ self.parameters.eB)

    def compute_public_key(self, private_key):
        pr = self.parameters
        KB = pr.PB + private_key * pr.QB
        phiB = pr.curve.isogeny(KB, algorithm="factored")
        return ( phiB.codomain(), phiB(pr.PA), phiB(pr.QA) )

    def compute_shared_secret(self, private_key, other_public_key):
        LB = other_public_key[1] + private_key * other_public_key[2]
        psiB = other_public_key[0].isogeny(LB, algorithm="factored")
        return psiB.codomain().j_invariant()

class SIDH_Parameters:
    def __init__(self, lA, lB, eA, eB, p, curve, f=1):
        self.lA = lA
        self.lB = lB
        self.eA = eA
        self.eB = eB
        self.p = p
        self.curve = curve

        # Choose points PA, QA, PB, QB
        self.PA, self.QA = ( (lB ** eB) * G  * f for G in curve.gens() )
        self.PB, self.QB = ( (lA ** eA) * G  * f for G in curve.gens() )

        # Verifications
        assert self.verify(), "SIDH parameters are not valid"

    def verify(self):
        curve = self.curve
        p = self.p
        lA, lB, eA, eB = self.lA, self.lB, self.eA, self.eB
        PA, QA, PB, QB = self.PA, self.QA, self.PB, self.QB
        print(f"{Back.LIGHTBLUE_EX}==== Verifying SIDH parameters [{self.__class__.__name__}] ==== {Style.RESET_ALL}")
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
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
            return False
        
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
        
        # points are of order lA^eA and lB^eB
        pa_correct_order = lA ** eA * PA == curve(0, 1, 0)
        qa_correct_order = lA ** eA * QA == curve(0, 1, 0)
        pb_correct_order = lB ** eB * PB == curve(0, 1, 0)
        qb_correct_order = lB ** eB * QB == curve(0, 1, 0)

        print(f"PA is of order lA^eA: {Back.LIGHTGREEN_EX if pa_correct_order else Back.RED}{pa_correct_order}{Style.RESET_ALL}")
        print(f"QA is of order lA^eA: {Back.LIGHTGREEN_EX if qa_correct_order else Back.RED}{qa_correct_order}{Style.RESET_ALL}")
        print(f"PB is of order lB^eB: {Back.LIGHTGREEN_EX if pb_correct_order else Back.RED}{pb_correct_order}{Style.RESET_ALL}")
        print(f"QB is of order lB^eB: {Back.LIGHTGREEN_EX if qb_correct_order else Back.RED}{qb_correct_order}{Style.RESET_ALL}")
        if not (pa_correct_order and qa_correct_order and pb_correct_order and qb_correct_order):
            print("Some points are not of the correct order")
            print(f"{Back.RED}==== SIDH parameters are not valid ==== {Style.RESET_ALL}")
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
    
        # Finished
        print(f"{Back.LIGHTBLUE_EX}==== SIDH parameters are valid ==== {Style.RESET_ALL}")
        return True


    def __str__(self):
        return f"SIDH parameters: \n lA: {self.lA} \n lB: {self.lB} \n eA: {self.eA} \n eB: {self.eB} \n p: {self.p} \n curve: {self.curve}"

class SIKEp182(SIDH_Parameters):
    def __init__(self):
        lA, eA, lB, eB = 2, 91, 3, 57
        p = (lA ** eA) * (lB ** eB) - 1
        F = FiniteField((p, 2), 'x', impl='pari_ffelt')
        curve = EllipticCurve(F, [1,0])
        super().__init__(lA, lB, eA, eB, p, curve)

class DeFeo(SIDH_Parameters):
    def __init__(self):
        lA, eA, lB, eB = 2, 63, 3, 41
        p = (lA ** eA) * (lB ** eB) * 11 - 1
        F = FiniteField((p, 2), 'x', impl='pari_ffelt')
        curve = EllipticCurve(F, [1,0])
        super().__init__(lA, lB, eA, eB, p, curve, 11)

class Baby_SIKE(SIDH_Parameters):
    def __init__(self):
        lA, eA, lB, eB = 2, 33, 3, 19
        p = (lA ** eA) * (lB ** eB) - 1
        F = FiniteField((p, 2), 'x', impl='pari_ffelt')
        curve = EllipticCurve(F, [1,0])
        super().__init__(lA, lB, eA, eB, p, curve)


class SIKEp434(SIDH_Parameters):
    def __init__(self):
        lA, eA, lB, eB = 2, 216, 3, 137
        p = (lA ** eA) * (lB ** eB) - 1
        F = FiniteField((p, 2), 'x', impl='pari_ffelt')
        curve = EllipticCurve(F, [1,0])
        super().__init__(lA, lB, eA, eB, p, curve)

class SIKEp751(SIDH_Parameters):
    def __init__(self):
        lA, eA, lB, eB = 2, 372, 3, 239
        p = (lA ** eA) * (lB ** eB) - 1
        F = FiniteField((p, 2), 'x', impl='pari_ffelt')
        curve = EllipticCurve(F, [1,0])
        super().__init__(lA, lB, eA, eB, p, curve)

available_curves ={
    "p182": SIKEp182,
    "DeFeo": DeFeo,
    "bSIKE": Baby_SIKE,
    "p434": SIKEp434,
    "p751": SIKEp751
}

def get_curve(curve_name):
    if curve_name not in available_curves:
        raise Exception(f"Curve {curve_name} not available")
    return available_curves[curve_name]()

def create_protocol(settings):
    partyA = SIDH_Party_A(settings)
    partyB = SIDH_Party_B(settings)
    return DH_Protocol(partyA, partyB)

