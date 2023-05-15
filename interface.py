# ==============================================================================
# Interface definition for DH type protocols
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

from sage.all import *
from colorama import Back, Style

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








