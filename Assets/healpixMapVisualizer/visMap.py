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
    opts = parser.parse_args(args)

    hpx = hp.read_map(opts.input)
    hp.mollview(hpx, coord="G", norm='log', min=1)
    hp.graticule()
    plt.show()


if __name__ == "__main__":
    main()


