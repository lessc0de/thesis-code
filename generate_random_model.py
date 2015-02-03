'''
Generate a random HMM model specification. Note that all entries in the
matrices will be greater than zero.
'''

import sys
import random
import argparse

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("--no_states", "-n", nargs=1, default=[7], type=int,
                        help="The number of states.")

    PARSER.add_argument("--alphabet", "-a", nargs=1, default=[4], type=int,
                        help="The size of the alphabet.")

    ARGS = PARSER.parse_args()

    NO_STATES = ARGS.no_states[0]
    ALPHABET_SIZE = ARGS.alphabet[0]

    PI = [random.random() for _ in range(NO_STATES)]
    PI = [x / sum(PI) for x in PI]

    A = [[random.random() for _ in range(NO_STATES)] for _ in range(NO_STATES)]
    for i in range(len(A)):
        A[i] = [x / sum(A[i]) for x in A[i]]

    B = [[random.random() for _ in range(ALPHABET_SIZE)] for _ in range(NO_STATES)]
    for i in range(len(B)):
        B[i] = [x / sum(B[i]) for x in B[i]]

    print("no_states")
    print(NO_STATES)
    print("alphabet_size")
    print(ALPHABET_SIZE)
    print("pi")
    for e in PI:
        print(e)
    print("A")
    for row in A:
        for e in row:
            sys.stdout.write(str(e) + " ")
        print()
    print("B")
    for row in B:
        for e in row:
            sys.stdout.write(str(e) + " ")
        print()
