"""
Support for manipulating Monce Carlo photon lists
"""

# standard python modules
import numpy as np
import struct

#SWD: Maybe should be rewritten as dictionary
class photons:
    """
    Class for storing photon list data
    """
    #Initialization from dictionary
    def __init__(self, phlist):
        self.npars = 10
        self.polarized = phlist['polarized']
        self.ntot = phlist['ntot']
        self.coord = phlist['coord']
        if (self.polarized):
            self.npars = self.npars + 2
        ncol = phlist['npars']
        self.nphot = phlist['length']
        if (ncol < self.npars):
            raise ValueError("Error creating photon: ncol {:d} < npars {:d}".format(ncol,self.npars))
        else:
            if (ncol > self.npars):
                self.nuser = ncol - self.npars
                self.user = np.zeros((self.nphot,self.nuser))
            else:
                self.nuser = 0

        # allocate arrays for each variable
        self.weight = phlist['list'][:,0]
        self.energy = phlist['list'][:,1]
        self.x0 = phlist['list'][:,5]
        self.x1 = phlist['list'][:,2]
        self.x2 = phlist['list'][:,3]
        self.x3 = phlist['list'][:,4]
        self.k0 = phlist['list'][:,9]
        self.k1 = phlist['list'][:,6]
        self.k2 = phlist['list'][:,7]
        self.k3 = phlist['list'][:,8]
        if (self.polarized):
            self.q = phlist['list'][:,10]
            self.u = phlist['list'][:,11]
            #self.v = phlist['list'][:,12]
        if (self.nuser > 0):
            for i in range(self.nuser):
                self.user[:,i] = phlist['list'][:,i+self.npars]

def read_list(filename):
    """
    Read unformated list output and return as a dictionary
    """

    try:
        # Read raw data
        with open(filename, 'rb') as data_file:
            raw_data = data_file.read()
        raw_data_ascii = raw_data.decode('ascii', 'replace')
    except:
        return None

    # Store in dictionary
    phlist = {}
    current_index = 0

    # Function for skipping though the file
    def skip_string(expected_string):
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index+expected_string_len] != expected_string:
            raise RuntimeError('File not formatted as expected')
        return current_index+expected_string_len

    current_index = skip_string("length=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['length'] = \
      list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("npars=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['npars'] = \
      list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("ntot=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['ntot'] = \
      list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("polarized=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['polarized'] = \
      bool(list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0])
    current_index = end_of_line_index + 1
    current_index = skip_string("coord=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['coord'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1

    # Read in data
    nelements = phlist['length'] * phlist['npars']
    format_string = '>' + 'd'*nelements
    begin_index = current_index
    end_index = begin_index + 8*nelements
    vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    phlist['list'] = vals.reshape((phlist['length'],phlist['npars']))

    return phlist

def write_list(filename,phlist):
    """
    Write photon list (dictionary) to file
    """
    outfile = open(filename, 'w')

    # Write header information
    outfile.write("length={:d}\n".format(phlist['length']))
    outfile.write("npars={:d}\n".format(phlist['npars']))
    outfile.write("ntot={:d}\n".format(phlist['ntot']))
    outfile.write("polarized={:d}\n".format(int(phlist['polarized'])))
    outfile.write("coord="+phlist['coord']+"\n")
    outfile.close()

    # Write list data
    outfile = open(filename, 'ab')
    nelements = phlist['length']*phlist['npars']
    phlist['list'] = phlist['list'].reshape(nelements)
    myfmt='>'+'d'*nelements
    bin=struct.pack(myfmt,*(phlist['list']))
    outfile.write(bin)
    outfile.close()
