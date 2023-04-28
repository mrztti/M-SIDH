# ==============================================================================
# Main runner entry point for the project
# Implemented using SAGE math library
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

import interface
import sage.all as sage


# ==============================================================================
# Attack on SIDH inspired from https://github.com/jack4818/Castryck-Decru-SageMath
# ==============================================================================

def attack_SIDH(public_params, target_secret):
    pass



if __name__ == "__main__":

    # ==============================================================================
    # TEST SIDH
    # ==============================================================================

    settings = interface.Baby_SIKE()
    partyA = interface.SIDH_Party_A(settings)
    partyB = interface.SIDH_Party_B(settings)

    protocol = interface.DH_Protocol(partyA, partyB)

    # Run protocol
    protocol.run()
