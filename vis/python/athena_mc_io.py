"""
Read in Athena++ monte carlo output data files.
"""

# standard python modules
import numpy as np
import struct

class photons:
    """
    Class for storing photon list data
    """
    #Initialization from dictionary
    def __init__(self, phlist):
        self.npars = 8
        self.polarized = phlist['polarized']
        self.relativistic = phlist['relativistic']
        self.ntot = phlist['ntot']
        self.coord = phlist['coord']
        if (self.polarized):
            self.npars = self.npars + 3
        if (self.relativistic):
            self.npars = self.npars + 2
        ncol = phlist['npars']
        self.nphot = phlist['length']
        if (ncol < self.npars):
            raise ValueError("Error creating photon: ncol {:d} < npars {:d}".format(ncol,self.npars))
        else:
            if (ncol < self.npars):
                self.nuser = ncol - self.npars
                self.user = np.zeros((self.nphot,self.nuser))
            else:
                self.nuser = 0

        
        # allocate arrays for each variable
        self.weight = phlist['list'][:,0]
        self.energy = phlist['list'][:,1]
        if (self.relativistic):
            self.x0 = phlist['list'][:,2]
            self.x1 = phlist['list'][:,3]
            self.x2 = phlist['list'][:,4]
            self.x3 = phlist['list'][:,5]
            self.k0 = phlist['list'][:,6]
            self.k1 = phlist['list'][:,7]
            self.k2 = phlist['list'][:,8]
            self.k3 = phlist['list'][:,9]
            if (self.polarized):
                self.q = phlist['list'][:,10]
                self.u = phlist['list'][:,11]
                self.v = phlist['list'][:,12]
        else:
            self.x1 = phlist['list'][:,2]
            self.x2 = phlist['list'][:,3]
            self.x3 = phlist['list'][:,4]
            self.k1 = phlist['list'][:,5]
            self.k2 = phlist['list'][:,6]
            self.k3 = phlist['list'][:,7]  
            if (self.polarized):
                self.q = phlist['list'][:,8]
                self.u = phlist['list'][:,9]
                self.v = phlist['list'][:,10]
        if (self.nuser > 0):
            for i in range(self.nuser):
                self.user[:,i] = phlist['list'][:,i+self.npars]

def readlist(filename):
    """
    Read unformated list output and return as a dictionary
    """
    # Read raw data
    with open(filename, 'rb') as data_file:
        raw_data = data_file.read()
    raw_data_ascii = raw_data.decode('ascii', 'replace')

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
    phlist['length'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("npars=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['npars'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("ntot=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['ntot'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("polarized=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['polarized'] = bool(map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0])
    current_index = end_of_line_index + 1
    current_index = skip_string("relativistic=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phlist['relativistic'] = bool(map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0])
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

def writelist(filename,phlist):
    """
    Write photon list (dictionary) to file
    """
    outfile = open(filename, 'w')
    
    # Write header information
    outfile.write("length={:d}\n".format(phlist['length']))
    outfile.write("npars={:d}\n".format(phlist['npars']))
    outfile.write("ntot={:d}\n".format(phlist['ntot']))
    outfile.write("polarized={:d}\n".format(int(phlist['polarized'])))
    outfile.write("relativistic={:d}\n".format(int(phlist['relativistic'])))
    outfile.write("coord="+phlist['coord']+"\n")

    # Write list data
    nelements = phlist['length']*phlist['npars']
    phlist['list'] = phlist['list'].reshape(nelements)
    myfmt='>'+'d'*nelements
    bin=struct.pack(myfmt,*(phlist['list']))
    outfile.write(bin)
    outfile.close()

def write_spectrum(filename,spectrum):
    """
    Writes spectrum to output file
    """

    # Open outfile
    outfile = open(filename, 'w')

    nx = spectrum['nx']
    nmu = spectrum['nmu']
    nphi = spectrum['nphi']

    # Write header information
    outfile.write("nx={:d}\n".format(nx))
    outfile.write("nmu={:d}\n".format(nmu))
    outfile.write("nphi={:d}\n".format(nphi))
    outfile.write("ntot={:d}\n".format(spectrum['ntot']))
    outfile.write("nintens={:d}\n".format(spectrum['nintens']))
    outfile.write("units="+spectrum['xaxis']+"\n")
    outfile.write("yerror="+spectrum['yerror']+"\n")
    # Write binfaces
    myfmt='>'+'d'*(nx+1)
    bin=struct.pack(myfmt,*(spectrum['xfaces']))
    outfile.write(bin)
    myfmt='>'+'d'*(nmu+1)
    bin=struct.pack(myfmt,*(spectrum['mufaces']))
    outfile.write(bin)
    myfmt='>'+'d'*(nphi+1)
    bin=struct.pack(myfmt,*(spectrum['phifaces']))
    outfile.write(bin)        
    # Write data
    nelements = (spectrum['nintens']*nx*nmu*nphi)
    myfmt='>'+'d'*nelements
    bin=struct.pack(myfmt,*(spectrum['intensity'].reshape(nelements)))
    outfile.write(bin)
    if (spectrum['yerror'] == 'true'):
        bin=struct.pack(myfmt,*(spectrum['errors'].reshape(nelements)))
        outfile.write(bin)
    outfile.close()

def read_spectrum(filename):
    """
    Read spectrum and return as a dictionary
    """

    # Read raw data
    with open(filename, 'rb') as data_file:
        raw_data = data_file.read()
    raw_data_ascii = raw_data.decode('ascii', 'replace')

    spectrum = {}
    current_index = 0
    
    # Function for skipping though the file
    def skip_string(expected_string):
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index+expected_string_len] != expected_string:
            raise RuntimeError('File not formatted as expected')
        return current_index+expected_string_len

    current_index = skip_string("nx=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nx'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nmu=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nmu'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nphi=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nphi'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("ntot=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['ntot'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nintens=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nintens'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("units=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['units'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("yerror=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['yerror'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1

    # Read in faces
    nx = spectrum['nx']
    format_string = '>' + 'd'*(nx+1)
    begin_index = current_index
    end_index = begin_index + 8*(nx+1)
    spectrum['xfaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    nmu = spectrum['nmu']
    format_string = '>' + 'd'*(nmu+1)
    begin_index = end_index
    end_index = begin_index + 8*(nmu+1)
    spectrum['mufaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    nphi = spectrum['nphi']
    format_string = '>' + 'd'*(nphi+1)
    begin_index = end_index
    end_index = begin_index + 8*(nphi+1)
    spectrum['phifaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))

    # Read intensities
    nintens = spectrum['nintens']
    nelements = nintens*nx*nmu*nphi
    format_string = '>' + 'd'*nelements
    begin_index = end_index
    end_index = begin_index + 8*nelements
    vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    spectrum['intensity'] = vals.reshape((nintens,nphi,nmu,nx))
    if (spectrum['yerror'] == 'true'):
        begin_index = end_index
        end_index = begin_index + 8*nelements
        vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
        spectrum['errors'] = vals.reshape((nintens,nphi,nmu,nx))
    return spectrum
