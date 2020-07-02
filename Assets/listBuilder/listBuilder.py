#!/usr/bin/python3

import os
import sys
from datetime import datetime
from argparse import ArgumentParser

start = datetime.now()

def getListOfFiles(opts):
    fileList = []
    if opts.subdirs:
        dirs = [ opts.input + "/" + tmpDir for tmpDir in os.listdir(opts.input) if os.path.isdir(opts.input + "/" + tmpDir) ]
        for tmpdir in dirs:
            if opts.verbose:
                print('Searching for files in {}'.format(tmpdir))
            tmpfiles = [ file for file in os.listdir(tmpdir) if os.path.isfile(tmpdir + "/" + file) and file.endswith('.root') ]
            if tmpfiles:
                for file in tmpfiles:
                    fileList.append(file)
    else:
        tmpfiles = [ file for file in os.listdir(opts.input) if os.path.isfile(opts.input + "/" + file) and file.endswith('.root') ]
        if tmpfiles:
            for file in tmpfiles:
                fileList.append(file)

    return fileList

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Data list builder")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='Input DATA WD')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='Output list (text file)')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')
    parser.add_argument("-s", "--subdirs", dest='subdirs', default=False,
                        action='store_true', help='search in subdirs')
    
    opts = parser.parse_args(args)

    # Get list of ROOT files
    dataList = getListOfFiles(opts)
    
    # Write files to list
    with open(opts.output, "w") as myList:
        for elm in dataList:
            myList.write(elm)


if __name__ == "__main__":
    main()
    print(datetime.now()-start)