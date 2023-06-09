# ==============================================================================
# Main runner entry point for the project
# Implemented using SAGE math library
#
# Author: Malo RANZETTI
# Date: Spring 2023
# ==============================================================================

import os
import msidh 
import sidh
import sage.all as sage
import time
import argparse
import numpy as np

parser = argparse.ArgumentParser(
                    prog='M-SIDH Demo Runner',
                    description='This program is a demo built using SAGE of proposed countermeasures to the SIDH scheme.',
                    epilog='Written by M.Ranzetti')


def test_SIDH(curve, n_rounds=10):
    
    # ==============================================================================
    # TEST SIDH
    # NIST-LEVEL 5 (AES-256 security) 21s for the runtime
    # ==============================================================================
    print("Testing SIDH protocol...")

    scheme = sidh.create_protocol(sidh.get_curve(curve))
    results = []
    for i in range(n_rounds):
        print(f"Round {i+1}/{n_rounds}")
        results.append(scheme.run())

    print(f"Average time: {sum([r[1]for r in results])/n_rounds * 1e-9}s")
    print(f"Failure count: {sum([1 for r in results if not r[0]])}")

    average_time = sum([r[1]for r in results])/n_rounds * 1e-9
    std = np.std([r[1]for r in results])
    failure_count = sum([1 for r in results if not r[0]])

    data = {
        'settings': curve,
        'average_time': average_time,
        'std': std,
        'failure_count': failure_count
    }

    return data

def test_MSIDH(filename, n_rounds=10):


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
    scheme = msidh.create_protocol_from_file(filename)

    results = []
    for i in range(n_rounds):
        
        print(f"Round {i+1}/{n_rounds}")
        results.append(scheme.run())


    print(f"Average time: {sum([r[1]for r in results])/n_rounds * 1e-9}s")
    print(f"Failure count: {sum([1 for r in results if not r[0]])}")

    average_time = sum([r[1]for r in results])/n_rounds * 1e-9
    std = np.std([r[1]for r in results])
    failure_count = sum([1 for r in results if not r[0]])

    data = {
        'settings': filename.split('AES-')[1].split('.')[0],
        'average_time': average_time,
        'std': std,
        'failure_count': failure_count
    }

    return data

def gen_MSIDH128():
    msidh.create_g128_protocol()

def create_msidh(lam):
    msidh.create_protocol(msidh.MSIDHpArbitrary, lam) 

def output_data(filename, data):
    '''
    Write the data given as an array into csv format
    data: dict of the form {name: [data list]}
    '''

    # check if file exists
    if not os.path.exists(filename):
        with open(filename, 'w') as f:
            f.write(','.join(data.keys()) + '\n')

    with open(filename, 'a') as f:
        # add a line
        f.write(','.join([str(v) for v in data.values()]) + '\n')
    print(f"Data written to {filename}")


if __name__ == "__main__":
    parser.add_argument('-t', '--test', type=str, choices=['sidh', 'msidh'], help='Test to run (sidh, msidh)')
    parser.add_argument('-c', '--curve', type=str, choices=list(sidh.available_curves.keys()) ,help='Curve to use for SIDH')
    parser.add_argument('-f', '--file', type=str, help='File to use for MSIDH paramters')
    parser.add_argument('-r', '--rounds', type=int, default=10, help='Number of rounds to run tests for')
    parser.add_argument('-g', '--gen', type=int, help='generate MSIDH parameters for a given security level')
    parser.add_argument('-g128', '--gen128', action='store_true', help='generate MSIDH-128 parameters')
    args = parser.parse_args()

    if args.gen:
        create_msidh(args.gen)
    elif args.gen128:
        gen_MSIDH128()
    elif args.test == 'sidh':
        if not args.curve:
            print("Please provide a curve to use for SIDH using -c")
            exit(1)
        data = test_SIDH(args.curve, args.rounds)
        output_data("sidh_results.csv", data)

    elif args.test == 'msidh':
        if not args.file:
            print("Please provide a file to use for MSIDH using -f")
            print("You can generate a file using -g <security level>")
            exit(1)
        data = test_MSIDH(args.file, args.rounds)
        output_data("msidh_results.csv", data)
    else:
        print("Invalid arguments, use -h for help")
        exit(1)

