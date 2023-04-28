# ==============================================================================
# SIDH: Supersingular Isogeny Diffie-Hellman
# Implemented using SAGE math library
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

from sage.all import *
from colorama import Fore, Back, Style

# ==============================================================================
# Helper functions
# ==============================================================================

def check_secrets(secret1, secret2):
    print("========================================")
    if secret1 == secret2:
        print(f"{Back.GREEN}{Style.BRIGHT}SUCCESSFUL SECRET SHARING{Style.RESET_ALL}")
        print(f"Shared secret: {secret1}")
        # print tzpe info
        print(f"Type: {type(secret1)}")
    else:
        print(f"{Back.RED}{Style.BRIGHT}FAILED SECRET SHARING{Style.RESET_ALL}")
        print(f"Secret 1: {secret1}")
        print(f"Secret 2: {secret2}")

    print("========================================")



# ==============================================================================

# A Party is a participant in the Diffie-Hellman protocol
class Party:
    def __init__(self, interface, name):
        self.interface = interface
        self.name = name

    def generate_private_key(self):
        self.private_key = self.interface.generate_private_key()

    def compute_public_key(self):
        self.public_key = self.interface.compute_public_key(self.private_key)

    def register_public_key(self, other_public_key):
        self.other_public_key = other_public_key

    def compute_shared_secret(self):
        self.shared_secret = self.interface.compute_shared_secret(self.private_key, self.other_public_key)

    def __str__(self):
        return f"Party ref.: {self.name}\nPrivate key: {self.private_key}\nPublic key: {self.public_key}\nShared: {self.shared_secret}"
    

# ==============================================================================
# Simulate a network
# ==============================================================================
class Pipe:
    def __init__(self, partyA, partyB):
        self.partyA = partyA
        self.partyB = partyB
        self.transmitted_messages_to_A = []
        self.transmitted_messages_to_B = []

    def transmit_A_to_B(self, message):
        self.transmitted_messages_to_B.append(message)
        self.partyB.register_public_key(message)

    def transmit_B_to_A(self, message):
        self.transmitted_messages_to_A.append(message)
        self.partyA.register_public_key(message)

    def get_total_trasmitted_bytes(self):
        # TODO: implement
        return len(self.transmitted_messages_to_A) + len(self.transmitted_messages_to_B)

# Protocol class for abstract Diffie-Hellman
class DH_interface:
    def __init__(self):
        raise NotImplementedError
    
    def get_public_parameters(self):
        raise NotImplementedError
    
    def print_public_parameters(self):
        raise NotImplementedError

    def generate_private_key(self):
        raise NotImplementedError

    def compute_public_key(self, private_key):
        raise NotImplementedError
    
    def compute_shared_secret(self, private_key, other_public_key):
        raise NotImplementedError
    
    def __str__(self):
        return f"{Style.DIM}DH-interface >> {self.__class__.__name__} \n {self.print_public_parameters()} {Style.RESET_ALL}"


class DH_Protocol:
    def __init__(self, interfaceA, interfaceB):
        self.interfaceA = interfaceA
        self.interfaceB = interfaceB

    def __str__(self):
        return f"DH_Protocol >> \nInterface 1 -->\n{self.interfaceA}\nInterface 2 --> \n{self.interfaceB}"

    def run(self):

        print("Running Diffie-Hellman protocol...")
        print(self)

        # Create parties
        alice = Party(self.interfaceA, "Alice")
        bob = Party(self.interfaceB, "Bob")

        # Create network
        network = Pipe(alice, bob)

        # Generate private keys
        alice.generate_private_key()
        bob.generate_private_key()

        # Compute public keys
        alice.compute_public_key()
        bob.compute_public_key()

        # Exchange public keys
        network.transmit_A_to_B(alice.public_key)
        network.transmit_B_to_A(bob.public_key)

        # Compute shared secrets
        alice.compute_shared_secret()
        bob.compute_shared_secret()

        check_secrets(alice.shared_secret, bob.shared_secret)

        return alice, bob

# ==============================================================================
# Regular diffie-hellman (using SAGE)
# ==============================================================================

class DH(DH_interface):
    def __init__(self, group = Zmod(7919)):
        self.group = group
        self.p = group.order()
        self.gen = group.random_element()

    def get_public_parameters(self):
        return self.group, self.p, self.gen
    
    def print_public_parameters(self):
        return f"Public parameters: \n Group: {self.group} \n p (order): {self.p} \n Generator: {self.gen}"

    def generate_private_key(self):
        return self.group.random_element()

    def compute_public_key(self, private_key):
        return self.gen ** private_key

    def compute_shared_secret(self, private_key, other_public_key):
        return other_public_key ** private_key
    

# ==============================================================================
# Supersingular Isogeny Diffie-Hellman (SIDH)
# ==============================================================================
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
    def __init__(self, lA, lB, eA, eB, p, curve):
        self.lA = lA
        self.lB = lB
        self.eA = eA
        self.eB = eB
        self.p = p
        self.curve = curve

        # Choose points PA, QA, PB, QB
        self.PA, self.QA = ( (lB ** eB) * G for G in curve.gens() )
        self.PB, self.QB = ( (lA ** eA) * G for G in curve.gens() )

        # Verifications
        assert self.verify(), "SIDH parameters are not valid"

    def verify(self):
        curve = self.curve
        p = self.p
        lA, lB, eA, eB = self.lA, self.lB, self.eA, self.eB
        PA, QA, PB, QB = self.PA, self.QA, self.PB, self.QB
        print(f"{Back.LIGHTBLUE_EX}==== Verifying SIDH parameters [{self.__class__.__name__}] ==== {Style.RESET_ALL}")
        supersingular = curve.is_supersingular()
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
        super().__init__(lA, lB, eA, eB, p, curve)

class Baby_SIKE(SIDH_Parameters):
    def __init__(self):
        lA, eA, lB, eB = 2, 33, 3, 19
        p = (lA ** eA) * (lB ** eB) - 1
        F = FiniteField((p, 2), 'x', impl='pari_ffelt')
        curve = EllipticCurve(F, [1,0])
        super().__init__(lA, lB, eA, eB, p, curve)








