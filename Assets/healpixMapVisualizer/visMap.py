#!/usr/bin/env python3

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib import cm

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Healpix map visualizer")

    parser.add_argument("-m", "--map", type=str, dest='input', help='Input map')
    parser.add_argument("-l", "--log", dest='log', default=False, action='store_true', help='log scale (min value set to 1)')

    opts = parser.parse_args(args)

    hpx = hp.read_map(opts.input)
    if opts.log:
        hp.mollview(hpx, coord="G", norm='log', min = 1)
    else:
        hp.mollview(hpx, coord="G", norm='None')
    
    hp.graticule()
    plt.show()


if __name__ == "__main__":
    main()


