# ==============================================================================
# Interface definition for DH type protocols
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

from sage.all import *
from colorama import Back, Style
import time

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
        print("========================================")
        return True
    else:
        print(f"{Back.RED}{Style.BRIGHT}FAILED SECRET SHARING{Style.RESET_ALL}")
        print(f"Secret 1: {secret1}")
        print(f"Secret 2: {secret2}")
        print("========================================")
        return False

# ==============================================================================

# A Party is a participant in the Diffie-Hellman protocol
class Party:
    def __init__(self, interface, name):
        self.interface = interface
        self.name = name

    def generate_private_key(self):
        self.private_key = self.interface.generate_private_key()

    def compute_public_key(self):
        try:
            self.private_key
        except AttributeError:
            raise ValueError("Cannot compute public key without private key")
        self.public_key = self.interface.compute_public_key(self.private_key)

    def register_public_key(self, other_public_key):
        self.other_public_key = other_public_key

    def compute_shared_secret(self):
        try:
            self.private_key
        except AttributeError:
            raise ValueError("Cannot compute shared secret without private key")
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

    def generate_pk(self):
        self.alice.generate_private_key()
        self.bob.generate_private_key()

    def compute_pubk(self):
        self.alice.compute_public_key()
        self.bob.compute_public_key()

    def compute_ss(self):
        self.alice.compute_shared_secret()
        self.bob.compute_shared_secret()

    def run(self):

        print(f"{Back.YELLOW}{Style.BRIGHT}--++-- STARTING PROTOCOL --++--{Style.RESET_ALL}")

        # Create parties
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- CREATING PARTIES --++--{Style.RESET_ALL}")
        self.alice = Party(self.interfaceA, "Alice")
        self.bob = Party(self.interfaceB, "Bob")
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- PARTIES CREATED --++--{Style.RESET_ALL}")

        # Create network
        network = Pipe(self.alice, self.bob)

        TIME = time.time_ns()

        # Generate private keys
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- GENERATING PRIVATE KEYS --++--{Style.RESET_ALL}")
        timer_start = time.time_ns()
        self.generate_pk()
        print(f"Elapsed time: {(time.time_ns() - timer_start) / 1e9} s")
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- PRIVATE KEYS GENERATED --++--{Style.RESET_ALL}")

        # Compute public keys
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- COMPUTING PUBLIC KEYS --++--{Style.RESET_ALL}")
        timer_start = time.time_ns()
        self.compute_pubk()
        print(f"Elapsed time: {(time.time_ns() - timer_start) / 1e9} s")
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- PUBLIC KEYS COMPUTED --++--{Style.RESET_ALL}")

        # Exchange public keys
        print(f"{Back.BLUE}{Style.BRIGHT}--++-- EXCHANGING PUBLIC KEYS --++--{Style.RESET_ALL}")
        network.transmit_A_to_B(self.alice.public_key)
        network.transmit_B_to_A(self.bob.public_key)


        # Compute shared secrets
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- COMPUTING SHARED SECRETS --++--{Style.RESET_ALL}")
        timer_start = time.time_ns()
        self.compute_ss()
        print(f"Elapsed time: {(time.time_ns() - timer_start) / 1e9} s")
        print(f"{Back.LIGHTBLACK_EX}{Style.BRIGHT}--++-- SHARED SECRETS COMPUTED --++--{Style.RESET_ALL}")


        total_time = time.time_ns() - TIME
        print(f"{Back.YELLOW}{Style.BRIGHT}--++-- PROTOCOL COMPLETED --++--{Style.RESET_ALL}")
        print(f"Total time: {total_time / 1e9} s")
        return check_secrets(self.alice.shared_secret, self.bob.shared_secret), total_time








