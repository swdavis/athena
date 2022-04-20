#! /usr/bin/env python

"""
Read in photon lists written by independent processes and join
into single file
"""

# python standard modules
import argparse
import numpy as np
from os import system

# Athena++ modules
import athena_mc as athenamc

# Main function
def main(**kwargs):

    nproc = kwargs['nproc']
    if kwargs['startzero']:
        pstart = 0
    else:
        pstart = 1

    start = kwargs['start']
    end = kwargs['end']

    basename = kwargs['basename']

    def list_match(list1,list2):
        match = True
        if (list1 is None):
            return False
        if (list2 is None):
            return False
        if (list1['npars'] != list2['npars']):
            match = False
        if (list1['polarized'] != list2['polarized']):
            match = False
        return match

    firstlist = True
    filelist = ""
    ntot = 0
    for i in range(pstart,nproc):
        for j in range(start,end+1):
            filename = basename+".proc{:d}".format(i)+".{:05d}".format(j)+".list"
            filelist += filename + " "
            if (firstlist):
                phlist = athenamc.read_list(filename)
                firstlist = False
                ntot += phlist['ntot']
            else:
                addlist = athenamc.read_list(filename)
                if (list_match(phlist,addlist)):
                    phlist['list'] = np.append(phlist['list'],addlist['list'])
                    phlist['length'] += addlist['length']
                    ntot += addlist['ntot']
                else:
                    if (addlist is None):
                        if (not kwargs["skip"]):
                            raise RuntimeError(filename+" not found, aborting.\n")
                    else:
                        raise RuntimeError("List headers do not match for "+filename+".\n")
            print("Reading: "+filename,ntot)
    phlist['ntot'] = ntot
    print("Final list contains {:d} photons out of {:d} initialized.\n"
          .format(phlist['length'],phlist['ntot']))
    athenamc.write_list(basename+".list",phlist)


    if kwargs['removeold']:
        # Remove source files after write
        system("rm "+filelist)
    else:
        # write filelist to file.  Allows manual remove once fidelity of the combined
        # list has been verified e.g. "rm $(cat basename.list.files)"
        listfile = open(basename+".list.files",'w')
        listfile.write(filelist)
        listfile.close()

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('basename',
        help='basename of list files')
    parser.add_argument('nproc',
        type=int,
        help='number of processes')
    parser.add_argument('start',
        type=int,
        help='starting output number')
    parser.add_argument('end',
        type=int,
        help='ending output number')
    parser.add_argument('--skip',
        action='store_true',
        help='include lists from first process')
    parser.add_argument('-s', '--startzero',
        action='store_true',
        help='include lists from first process')
    parser.add_argument('-rm', '--removeold',
        action='store_true',
        help='delete origin files')

    args = parser.parse_args()
    main(**vars(args))
