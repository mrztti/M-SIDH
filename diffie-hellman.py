# ==============================================================================
# Regular Diffie-Hellman key exchange
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

from sage.all import *
from interface import DH_interface

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