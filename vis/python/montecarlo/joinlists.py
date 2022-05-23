#! /usr/bin/env python

"""
Read in photon lists written by independent processes and join
into single file
"""

#SWD modify writing to file in the case of larger of memory cap

# python standard modules
import argparse
import numpy as np
from os import system

# Athena++ modules
import athena_mc as athenamc

# Check if headers of lists match
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

def join(filelist,outfile,skip):

    # First Read through and check headers match, get ntot
    firstlist = True
    ntot = 0
    length = 0
    for infile in filelist:
        if (firstlist):
            phlist = athenamc.read_list(infile,data=False)
            firstlist = False
            ntot += phlist['ntot']
            length += phlist['ntot']
        else:
            addlist = athenamc.read_list(infile,data=False)
            if (list_match(phlist,addlist)):
                ntot += addlist['ntot']
                length += addlist['length']
            else:
                if (addlist is None):
                    if (not skip):
                        raise RuntimeError(infile+" not found, aborting.\n")
                else:
                    raise RuntimeError("List headers do not match for "+infile+".\n")
    print("Final list contains {:d} photons out of {:d} initialized.\n"
          .format(length,ntot))

    # If headers all match, read data and write output file as
    firstlist = True
    for infile in filelist:
        print("Reading: "+infile)
        if (firstlist):
            phlist = athenamc.read_list(infile,data=True)
            phlist['ntot'] = ntot
            firstlist = False
            athenamc.write_list(outfile,phlist,header=True,length=length)
        else:
            phlist = athenamc.read_list(infile,data=True)
            athenamc.write_list(outfile,phlist,header=False)

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
    outfile = basename+'.list'

    fileliststring = ""
    filelist = []
    for j in range(start,end+1):
        for i in range(pstart,nproc):
            filename = basename+".proc{:d}".format(i)+".{:05d}".format(j)+".list"
            filelist.append(filename)
            fileliststring += filename + " "

    join(filelist,outfile,kwargs['skip'])

    if kwargs['removeold']:
        # Remove source files after write
        system("rm "+fileliststring)
    else:
        # write filelist to file.  Allows manual remove once fidelity of the combined
        # list has been verified e.g. "rm $(cat basename.list.files)"
        listfile = open(basename+".list.files",'w')
        listfile.write(fileliststring)
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
