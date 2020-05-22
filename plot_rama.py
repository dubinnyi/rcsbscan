#!/usr/bin/python3 -u

import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    arg_parser = argparse.ArgumentParser(
        usage='Plot phi,psi on ramachandran from two numpy arrays')
    arg_parser.add_argument('phi', type=str, help='numpy ndarray with phi angles')
    arg_parser.add_argument('psi', type=str, help='numpy ndarray with psi angles')
    args = arg_parser.parse_args()

    phi = np.load(args.phi)
    psi = np.load(args.psi)
    print("phi shape {}".format(phi.shape))
    print("psi shape {}".format(psi.shape))
    ndata = min(phi.shape[0], psi.shape[0])
    nres = min(phi.shape[1], psi.shape[1])
    for res in range(nres):
        plt.plot(phi[:,res], psi[:,res], 'b.', alpha=0.2)
        plt.xlabel('phi')
        plt.ylabel('psi')
        plt.title("Residue {}".format(res+1))
        plt.show()


if __name__ == "__main__":
    main()