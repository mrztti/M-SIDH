# ==============================================================================
# Main runner entry point for the project
# Implemented using SAGE math library
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

import msidh 
import sidh
import sage.all as sage


# ==============================================================================
# Attack on SIDH inspired from https://github.com/jack4818/Castryck-Decru-SageMath
# ==============================================================================

def attack_SIDH(public_params, target_secret):
    pass



if __name__ == "__main__":

    """ # ==============================================================================
    # TEST SIDH
    # ==============================================================================
    print("Testing SIDH protocol...")
    ex = sidh.create_protocol(sidh.get_curve("DeFeo"))

    # Run protocol
    ex.run() """

    # ==============================================================================
    # TEST MSIDH
    # ==============================================================================
    print("Testing MSIDH protocol...")
    settings = msidh.MSIDHpBaby()
    ex = msidh.create_protocol(settings).run()
