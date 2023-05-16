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
    sidh.create_protocol(sidh.get_curve("DeFeo")).run()
    """

    # ==============================================================================
    # TEST MSIDH
    # Current maximum tested: t = 90 
    # Settings generation: 91.5s
    # Protocol execution: 29.2s
    # Currently the biggest bottlenecks are: 
    # - prime verification in EllipticCurve
    # - computing the generators of the curve
    # 
    # ==============================================================================
    print("Testing MSIDH protocol...")
    settings = msidh.MSIDHp32
    msidh.create_protocol(settings).run()
