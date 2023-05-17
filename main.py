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


if __name__ == "__main__":

    # ==============================================================================
    # TEST SIDH
    # NIST-LEVEL 5 (AES-256 security) 21s for the runtime
    # ==============================================================================
    """ print("Testing SIDH protocol...")
    sidh.create_protocol(sidh.get_curve("p751")).run() """

    # ==============================================================================
    # TEST MSIDH
    # Current maximum tested: t = 90 // 100 // 200
    # GOAL -> t = 572 for AES-128 security
    # Settings generation: 32.8s // 294.6s // 194.4s
    # Protocol execution: 5.3s // 35.0s // 320.9s
    # Currently the biggest bottlenecks are: 
    # - prime verification in EllipticCurve (OVERRIDEN IN SAGE SOURCE CODE)
    # - computing the generators of the curve (=> there might be a way to optimize this)
    # 
    # ==============================================================================
    print("Testing MSIDH protocol...")
    settings = msidh.MSIDHpArbitrary
    msidh.create_protocol(settings, 64).run()
