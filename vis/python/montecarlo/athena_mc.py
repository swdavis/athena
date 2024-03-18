"""
Support for manipulating and plotting Monte Carlo outputs
"""

# standard python modules
import numpy as np
import struct
import gc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
import gc
import math

#SWD: Maybe photons be rewritten simply as dictionary
#SWD: Add error control
class photons:
    """
    Class for storing photon list data
    """
    #Initialization from dictionary
    def __init__(self, phlist):
        self.npars = 10
        self.dt = phlist['dt']
        self.polarized = phlist['polarized']
        self.ntot = phlist['ntot']
        self.coord = phlist['coord']
        if (self.polarized):
            self.npars = self.npars + 2
        ncol = phlist['npars']
        self.nphot = phlist['length']
        if (ncol < self.npars):
            raise ValueError("Error creating photon: ncol {:d} < npars {:d}"
                             .format(ncol,self.npars))
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

def read_list(filename,data=True,header=True):
    """
    Read unformated list output and return as a dictionary
    """

    mxh_ = 1000
    mxl_ = 1000000

    try:
        # Read raw data
        with open(filename, 'rb') as data_file:
            raw_data = data_file.read()
        raw_data_ascii = raw_data[0:mxh_].decode('ascii', 'replace')
    except:
        print("Could not open "+filename+" for reading.")
        return None

    # Store in dictionary
    phlist = {}
    current_index = 0

    # Function for skipping though the file
    def skip_string(expected_string):
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index+expected_string_len] != \
           expected_string:
            raise RuntimeError('File not formatted as expected')
        return current_index+expected_string_len

    if (header):
        try:
            current_index = skip_string("dt=")
            end_of_line_index = current_index + 1
            while raw_data_ascii[end_of_line_index] != '\n':
                end_of_line_index += 1
            phlist['dt'] = list(map(float,
                           raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
            current_index = end_of_line_index + 1
        except:
            print("List file contains no dt entry. Setting to 1.")
            phlist['dt'] = 1.

        current_index = skip_string("length=")
        end_of_line_index = current_index + 1
        while raw_data_ascii[end_of_line_index] != '\n':
            end_of_line_index += 1
        phlist['length'] = list(map(int,
                           raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
        current_index = end_of_line_index + 1

        current_index = skip_string("npars=")
        end_of_line_index = current_index + 1
        while raw_data_ascii[end_of_line_index] != '\n':
            end_of_line_index += 1
        phlist['npars'] = list(map(int,
            raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
        current_index = end_of_line_index + 1

        current_index = skip_string("ntot=")
        end_of_line_index = current_index + 1
        while raw_data_ascii[end_of_line_index] != '\n':
            end_of_line_index += 1
        phlist['ntot'] = list(map(int,
                         raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
        current_index = end_of_line_index + 1

        current_index = skip_string("polarized=")
        end_of_line_index = current_index + 1
        while raw_data_ascii[end_of_line_index] != '\n':
            end_of_line_index += 1
        phlist['polarized'] = bool(list(map(int,
            raw_data_ascii[current_index:end_of_line_index].split(' ')))[0])
        current_index = end_of_line_index + 1

        current_index = skip_string("coord=")
        end_of_line_index = current_index + 1
        while raw_data_ascii[end_of_line_index] != '\n':
            end_of_line_index += 1
        phlist['coord'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
        current_index = end_of_line_index + 1

    old = False
    if (data):
        # Read in data
        npars = phlist['npars']
        length = phlist['length']
        nelements = length * npars

        begin_index = current_index
        sizeloop = mxl_*npars
        nloop = nelements // sizeloop
        sizelast = nelements % sizeloop
        # SWD: temporary to read broken list
        #for i in range(nloop):
        ti = time.time()
        for i in range(nloop+1):
            vals = ()
            if (i == nloop):
                size = sizelast
            else:
                size = sizeloop
            end_index = begin_index + 8*size
            format_string = '>' + 'd'*size
            if (i > 0):
                print(i,"/",nloop)
                tf = time.time()
                print(tf-ti)
                ti = tf
            # Check if end_index is larger than size of raw_data
            # if so, resize list
            lrd = len(raw_data)
            if lrd < end_index:
                deficit = math.ceil((end_index - lrd) / (8*npars))
                print("Warning raw_data is smaller than expected by at least"
                      " {:d} samples.".format(deficit))
                end_index -= deficit*npars*8
                length -= deficit
                phlist['length'] = length
                format_string = '>' + 'd'*(size-deficit*npars)
            vals = struct.unpack(format_string, raw_data[begin_index:end_index])
            begin_index = end_index
            if (i == 0):
                phlist['list'] = np.array(vals)
            else:
                phlist['list'] = np.append(phlist['list'],np.array(vals))
            del vals
            gc.collect()
        phlist['list'] = phlist['list'].reshape((length,npars))
        # SWD: temporary to read broken list
        #phlist['length'] = nloop*mxl_
        #phlist['list'] = phlist['list'].reshape((nloop*mxl_,npars))
    #return phlist, current_index
    return phlist

def write_list(filename,phlist,header=True,length=None):
    """
    Write photon list (dictionary) to file
    """
    if (header):
        outfile = open(filename, 'w')
        # Write header information
        outfile.write("dt={:.8e}\n".format(phlist['dt']))
        if (length is None):
            outfile.write("length={:d}\n".format(phlist['length']))
        else:
            outfile.write("length={:d}\n".format(length))
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

def get_luminosity_list(phlist):
    """
    Read in list file and compute luminoisty
    """
    phots = photons(phlist)
    lumin = np.sum(phots.weight*phots.energy)/phots.dt

    return lumin

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
    outfile.write("dt={:.8e}\n".format(spectrum['dt']))
    outfile.write("nx={:d}\n".format(nx))
    outfile.write("nmu={:d}\n".format(nmu))
    outfile.write("nphi={:d}\n".format(nphi))
    outfile.write("ntot={:d}\n".format(spectrum['ntot']))
    outfile.write("nintens={:d}\n".format(spectrum['nintens']))
    outfile.write("units="+spectrum['xaxis']+"\n")
    outfile.write("polarized="+spectrum['polarized']+"\n")
    outfile.write("yerror="+spectrum['yerror']+"\n")
    outfile.close()

    # Write binfaces
    outfile = open(filename, 'ab')
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
            print(raw_data_ascii[current_index:current_index+expected_string_len])
            raise RuntimeError('File not formatted as expected')
        return current_index+expected_string_len

    try:
        current_index = skip_string("dt=")
        end_of_line_index = current_index + 1
        while raw_data_ascii[end_of_line_index] != '\n':
            end_of_line_index += 1
        spectrum['dt'] = list(map(float,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
        current_index = end_of_line_index + 1
    except:
        print("Spectrum file contains no dt entry. Setting to 1.")
        spectrum['dt'] = 1.

    current_index = skip_string("nx=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nx'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nmu=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nmu'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nphi=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nphi'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("ntot=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['ntot'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nintens=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nintens'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("units=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['units'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("polarized=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['polarized'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
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
    spectrum['xfaces'] = np.array(struct.unpack(format_string,
                                                raw_data[begin_index:end_index]))
    nmu = spectrum['nmu']
    format_string = '>' + 'd'*(nmu+1)
    begin_index = end_index
    end_index = begin_index + 8*(nmu+1)
    spectrum['mufaces'] = np.array(struct.unpack(format_string,
                                                 raw_data[begin_index:end_index]))
    nphi = spectrum['nphi']
    format_string = '>' + 'd'*(nphi+1)
    begin_index = end_index
    end_index = begin_index + 8*(nphi+1)
    spectrum['phifaces'] = np.array(struct.unpack(format_string,
                                                  raw_data[begin_index:end_index]))

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

# Check if headers match
def header_match(dict1,dict2,dict_type):

    match = True
    if (dict1 is None):
        return False
    if (dict2 is None):
        return False
    if (dict_type == 'list'):
        if (dict1['npars'] != dict2['npars']):
            match = False
        elif (dict1['polarized'] != dict2['polarized']):
            match = False
    elif (dict_type == 'spec'):
        if (dict1['nx'] != dict2['nx']):
            match = False
        elif (dict1['nmu'] != dict2['nmu']):
            match = False
        elif (dict1['nphi'] != dict2['nphi']):
            match = False
        elif (dict1['nintens'] != dict2['nintens']):
            match = False
        elif (dict1['units'] != dict2['units']):
            match = False
        elif (dict1['yerror'] != dict2['yerror']):
            match = False
    else:
        print("file type: "+dict_type+" not supported. Returning false.")
        match = False

    return match

def add_spectra(spec1,spec2,method='statistical'):
    """
    Add two spectra to create single spectrum
    """
    if (not header_match(spec1,spec2,'spec')):
        raise RuntimeError('[add_specta]: headers do not match')

    # copy spectra to new dictionaries to avoid modification of originals
    spec1c = spec1.copy()
    spec2c = spec2.copy()

    if (method == 'statistical'):
        if (spec1c['dt'] != spec2c['dt']):
            raise RuntimeError('[add_specta]: statistical averaging requested' \
                               'but spectra have different integration times')
        # undo normalization by photon count
        spec1c['intensity'] *= spec1c['ntot']
        spec2c['intensity'] *= spec2c['ntot']
        if (spec1c['yerror']):
            spec1c['errors'] *= spec1c['ntot']
            spec2c['errors'] *= spec2c['ntot']
    #elif (method == 'temporal'):
    #    if (spec1c['ntot'] != spec2c['ntot']):
    #        raise RuntimeError('[add_specta]: temporal averaging requested' \
    #                           'but spectra have different total counts')
    #    # undo normalization by integration time
    #    spec1c['intensity'] *= spec1c['dt']
    #    spec2c['intensity'] *= spec2c['dt']
    #    if (spec1c['yerror']):
    #        spec1c['errors'] *= (spec1c['dt'])**2
    #        spec1c['errors'] *= (spec2c['dt'])**2

    # intialize spec_out as copy for simplicity
    spec_out = spec1c.copy()

    # sum unormalized intensity and error arrays
    spec_out['intensity'] = np.add(spec1c['intensity'],spec2c['intensity'])
    if (spec_out['yerror']):
        spec_out['errors'] = np.add(spec1c['errors'],spec2c['errors'])

    spec_out['ntot'] = spec1c['ntot']+ spec2['ntot']

    if (method == 'statistical'):
        # renormalization by photon number
        spec_out['intensity'] /= spec_out['ntot']
        if (spec_out['yerror']):
            spec_out['errors'] /= spec_out['ntot']


    return spec_out

# Retrun x for desired units
def convert_xaxis(newunit,spectrum):

    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10

    baseunit = spectrum['units']
    xfaces = spectrum['xfaces']
    if (baseunit == 'kev'):
        nu = xfaces*1000.*everg/h
    if (baseunit == 'ev'):
        nu = xfaces*everg/h
    if (baseunit == 'nu'):
        nu = xfaces
    if (baseunit == 'lambda'):
        nu = c/1.e8/xfaces
    if (newunit == 'kev'):
        spectrum['xfaces'] = nu*h/(everg*1000.)
    if (newunit == 'ev'):
        spectrum['xfaces'] = nu*h/everg
    if (newunit == 'nu'):
        spectrum['xfaces'] = nu
    if (newunit == 'lambda'):
        spectrum['xfaces'] = c/nu*1.e8
    spectrum['units'] = newunit

# Return nu
def get_frequency(xunit,xfaces):

    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10

    if (xunit == 'kev'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])*1000.*everg/h
    if (xunit == 'ev'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])*everg/h
    if (xunit == 'nu'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])
    if (xunit == 'lambda'):
        nu = 0.5*(1./xfaces[:-1]+1./xfaces[1:])*c/1.e8
    return nu

def compute_nulnu_error(intensity,nu,errors=None):
    """
    Compute nu*L_nu and error
    """

    y = intensity[0,:]*nu
    if errors is not None:
        err = errors[0,:]*nu
        return y, err
    else:
        return y, None

def compute_lnu_error(intensity,errors=None):
    """
    Compute L_nu and error
    """
    y = intensity[0,:]
    if errors is not None:
        yerr = errors[0,:]
        return y, yerr
    else:
        return y, None

def compute_counts_error(intensity,nu,errors=None):
    """
    Compute L_nu and error
    """
    h = 6.62607015e-27
    y = intensity[0,:]/(h*nu)
    if errors is not None:
        yerr = errors[0,:]/(h*nu)
        return y, yerr
    else:
        return y, None

def compute_pol_frac_error(intensity,errors=None):
    """
    Compute fractional polarization and error (%)
    """
    i = intensity[0,:]
    q = intensity[1,:]
    u = intensity[2,:]
    frac = np.sqrt(q**2+u**2)/i
    if errors is not None:
        ei = errors[0,:]
        eq = errors[1,:]
        eu = errors[2,:]
        err = np.sqrt(((q*eq)**2+(u*eu)**2)/(q**2+u**2) + (q**2+u**2)*(ei/i)**2)/i
        return frac*100, err*100.
    else:
        return frac*100, None

def compute_pol_angle_error(intensity,errors=None):
    """
    Compute polarization angle and error if requested (degrees)
    """
    q = intensity[1,:]
    u = intensity[2,:]
    angle = 90./np.pi*np.arctan2(u,q)
    if errors is not None:
        eq = errors[1,:]
        eu = errors[2,:]
        err = 90./np.pi*np.sqrt(((u*eu)**2+(q*eq)**2))/(q**2+u**2)
        return angle, err
    else:
        return angle, None

def compute_q_error(intensity,errors=None):
    """
    Compute q=-Q/I and error if requested
    """
    i = intensity[0,:]
    q = intensity[1,:]
    frac = -q/i
    if errors is not None:
        ei = errors[0,:]
        eq = errors[1,:]
        err = np.sqrt(eq**2 + (q**2)*(ei/i)**2)/i
        return frac, err
    else:
        return frac, None

def compute_u_error(intensity,errors=None):
    """
    Compute u=U/I and error if requested
    """
    i = intensity[0,:]
    u = intensity[2,:]
    frac = u/i
    if errors is not None:
        ei = errors[0,:]
        eu = errors[1,:]
        err = np.sqrt(eu**2 + (u**2)*(ei/i)**2)/i
        return frac, err
    else:
        return frac, None

def compute_flux_frac_error(intensity,mufaces,errors=None):
    """
    Compute I/F and error if requested
    """
    i = intensity[0,:]
    mu = 0.5*(mufaces[1:]+mufaces[:-1])
    dmu = mufaces[1:]-mufaces[:-1]
    flux = np.sum(mu*dmu*i)
    frac = i/flux/2.
    if errors is not None:
        ei = errors[0,:]
        err = np.sqrt(ei**2 + np.sum((dmu*mu*ei)**2)*i**2/flux**2)/flux/2.
        return frac, err
    else:
        return frac, None

def polarization_requested(yunit):

    if (yunit == 'polfrac'):
        return True
    elif (yunit == 'polangle'):
        return True
    elif (yunit == 'q'):
        return True
    elif (yunit == 'u'):
        return True
    else:
        return False

def plot_frequency(spectrum,imu,iphi='ave',xunit='kev',yunit='nulnu',
                   plterr=True,nu=None):
    """
    Generate plot versus frequency (or equivalent).
    """

    # Set up x axis as bin midpoints
    xfaces = spectrum['xfaces']
    x = 0.5*(xfaces[1:]+xfaces[:-1])
    if (nu is None):
        nu = get_frequency(spectrum['units'],xfaces)

    # Initialize x labels
    xlabel = None
    if (xunit == 'kev'):
        xlabel = r"$E {\rm (keV)}$"
    if (xunit == 'ev'):
        xlabel = r"$E {\rm (eV)}$"
    if (xunit == 'nu'):
        xlabel = r"$\nu {\rm (Hz)}$"
    if (xunit == 'lambda'):
        xlabel = r"$\lambda {\rm (\AA)$"

    # Check if error requested and stored
    if plterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            plterr = False

    # Check whether spectrum has required polarization data
    if (polarization_requested(yunit) and (spectrum['polarized'] != 'true')):
        print("Error: polarization output "+yunit+" requested for unpolarized spectrum.")
        return None

    # Compute intensity spectrum
    intensity = spectrum['intensity']
    if plterr:
        errors = spectrum['errors']
    else:
        errors = None

    # Selection for azimuthal angle
    if ((iphi == 'ave') or (iphi == 'sum')):
        norm = 1./float(spectrum['nphi'])
        if iphi == 'sum':
            norm *= 2.*np.pi
        intensity = np.sum(intensity,axis=1)*norm
        if plterr:
            errors = np.sqrt(np.sum((errors)**2,axis=1))*norm
    else:
        iphi = int(iphi)
        intensity = intensity[:,iphi,:,:]
        if plterr:
            errors = errors[:,iphi,:,:]

    # Selection for polar angle
    if imu == 'sum':
        nmu = spectrum['nmu']
        mumid = 0.5*(spectrum['mufaces'][1:]+spectrum['mufaces'][:-1])
        intensity = np.tensordot(mumid,intensity,axes=[0,1])/nmu
        if plterr:
            errors = np.sqrt(np.tensordot((mumid)**2,(errors)**2,axes=[0,1]))/nmu
    else:
        imu = int(imu)
        intensity = intensity[:,imu,:]
        if plterr:
            errors = errors[:,imu,:]

    # Set y, yerr, and ylabel according to input units
    yerr = None
    ylabel = None
    if (yunit == 'nulnu'):
        ylabel = r"$\nu L_\nu {\rm (erg/s)}$"
        y, yerr = compute_nulnu_error(intensity,nu,errors)
    elif (yunit == 'lnu'):
        ylabel = r"$L_\nu {\rm (erg/s/Hz)}$"
        y, yerr = compute_lnu_error(intensity,errors)
    elif (yunit == 'counts'):
        ylabel = r"$N_\nu {\rm (counts/s/Hz)}$"
        y, yerr = compute_counts_error(intensity,nu,errors)
    elif (yunit == 'polfrac'):
        ylabel = r"$\rm Pol.\; Fraction \; (\%)$"
        y, yerr = compute_pol_frac_error(intensity,errors)
    elif (yunit == 'polangle'):
        ylabel = r"$\rm Pol.\; Angle$"
        y, yerr = compute_pol_angle_error(intensity,errors)
    elif (yunit == 'q'):
        ylabel = r"$Q_\nu/I_\nu$"
        y, yerr = compute_q_error(intensity,errors)
    elif (yunit == 'u'):
        ylabel = r"$U_\nu/I_\nu$"
        y, yerr = compute_u_error(intensity,errors)
    else:
        print("Error: yunit ("+yunit+") not specified correctly")
        return None

    # Return x and y variables, their labels, and possible error on y
    return x,y,yerr,xlabel,ylabel

def plot_theta(spectrum,ix,iphi='ave',xunit='mu',yunit='lnu',
               plterr=True,nu=None,verbose=False):
    """
    Generate plot versus polar angle (theta)
    """

    # Set up x axis as bin midpoints
    xfaces = spectrum['mufaces']
    x = 0.5*(xfaces[1:]+xfaces[:-1])

    # Initialize x labels
    xlabel = None
    if (xunit == 'mu'):
        xlabel = r"$\cos \theta$"

    # Check if error requested and stored
    if plterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            plterr = False

    # Check whether spectrum has required polarization data
    if (polarization_requested(yunit) and (spectrum['polarized'] != 'true')):
        print("Error: polarization output "+yunit+" requested for unpolarized spectrum.")
        return None

    intensity = spectrum['intensity']
    if plterr:
        errors = spectrum['errors']
    else:
        errors = None

    # Selection for azimuthal angle
    if ((iphi == 'ave') or (iphi == 'sum')):
        norm = 1./float(spectrum['nphi'])
        if iphi == 'sum':
            norm *= 2.*np.pi
        intensity = np.sum(intensity,axis=1)*norm
        if plterr:
            errors = np.sqrt(np.sum((errors)**2,axis=1))*norm
    else:
        iphi = int(iphi)
        intensity = intensity[:,iphi,:,:]
        if plterr:
            errors = errors[:,iphi,:,:]

    # Selection for frequency
    if ix == 'sum':
        nx = spectrum['nx']
        dx = (spectrum['xfaces'][1:]-spectrum['xfaces'][:-1]).reshape(nx)
        #intensity = np.dot(dx,intensity,axis=1)
        intensity = np.tensordot(dx,intensity,axes=[0,2])
        if plterr:
            errors = np.sqrt(np.tensordot((dx)**2,(errors)**2,axes=[0,2]))
    else:
        ix = int(ix)
        intensity = intensity[:,:,ix]
        if plterr:
            errors = errors[:,:,ix]

    # Set y, yerr, and ylabel according to input units
    yerr = None
    ylabel = None
    if (yunit == 'nulnu'):
        ylabel = r"$\nu L_\nu {\rm (erg/s)}$"
        y, yerr = compute_nulnu_error(intensity,nu,errors)
    elif (yunit == 'lnu'):
        ylabel = r"$L_\nu {\rm (erg/s/Hz)}$"
        y, yerr = compute_lnu_error(intensity,errors)
    elif (yunit == 'counts'):
        ylabel = r"$N_\nu {\rm (counts/s/Hz)}$"
        y, yerr = compute_counts_error(intensity,nu,errors)
    elif (yunit == 'polfrac'):
        ylabel = r"$\rm Pol.\; Fraction \; (\%)$"
        y, yerr = compute_pol_frac_error(intensity,errors)
    elif (yunit == 'polangle'):
        ylabel = r"$\rm Pol.\; Angle$"
        y, yerr = compute_pol_angle_error(intensity,errors)
    elif (yunit == 'q'):
        ylabel = r"$Q_\nu/I_\nu$"
        y, yerr = compute_q_error(intensity,errors)
    elif (yunit == 'u'):
        ylabel = r"$U_\nu/I_\nu$"
        y, yerr = compute_u_error(intensity,errors)
    elif (yunit == 'fluxfrac'):
        ylabel = r"$I_\nu/F_\nu$"
        y, yerr = compute_flux_frac_error(intensity,xfaces,errors)

    if (verbose):
        print(yunit)
        if (isinstance(ix,int)):
            print("x: ",spectrum['xfaces'][ix],' ',spectrum['units'])
        else:
            print("x: ",ix)
        if (isinstance(iphi,int)):
            print("phi: ",spectrum['phifaces'][iphi])
        else:
            print("phi: ",iphi)

    # Return x and y variables, there labels, and possible error on y
    return x,y,yerr,xlabel,ylabel

def plot_phi(spectrum,ix,imu='sum',xunit='phi',yunit='lnu',
               plterr=True,nu=None):
    """
    Generate plot versus azimuthal angle (phi)
    """

    # Set up x axis as bin midpoints
    xfaces = spectrum['phifaces']/(2.*np.pi)
    x = 0.5*(xfaces[1:]+xfaces[:-1])

    # Initialize x labels
    xlabel = None
    if (xunit == 'phi'):
        xlabel = r"$\phi/(2\pi)$"

    # Check if error requested and stored
    if plterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            plterr = False

    # Check whether spectrum has required polarization data
    if (polarization_requested(yunit) and (spectrum['polarized'] != 'true')):
        print("Error: polarization output "+yunit+" requested for unpolarized spectrum.")
        return None

    intensity = spectrum['intensity']
    if plterr:
        errors = spectrum['errors']
    else:
        errors = None

    # Selection for polar angle
    if (imu == 'sum'):
        nmu = spectrum['nmu']
        mumid = 0.5*(spectrum['mufaces'][1:]+spectrum['mufaces'][:-1])
        intensity = np.tensordot(mumid,intensity,axes=[0,2])/nmu
        if plterr:
            errors = np.sqrt(np.tensordot((mumid)**2,(errors)**2,axes=[0,2]))/nmu
    else:
        imu = int(imu)
        intensity = intensity[:,:,imu,:]
        if plterr:
            errors = errors[:,:,imu,:]

    # Selection for frequency
    if ix == 'sum':
        nx = spectrum['nx']
        dx = (spectrum['xfaces'][1:]-spectrum['xfaces'][:-1]).reshape(nx)
        #intensity = np.dot(dx,intensity,axis=1)
        intensity = np.tensordot(dx,intensity,axes=[0,2])
        if plterr:
            errors = np.sqrt(np.tensordot((dx)**2,(errors)**2,axes=[0,2]))
    else:
        ix = int(ix)
        intensity = intensity[:,:,ix]
        if plterr:
            errors = errors[:,:,ix]

    # Set y, yerr, and ylabel according to input units
    yerr = None
    ylabel = None
    if (yunit == 'nulnu'):
        ylabel = r"$\nu L_\nu {\rm (erg/s)}$"
        y, yerr = compute_nulnu_error(intensity,nu,errors)
    elif (yunit == 'lnu'):
        ylabel = r"$L_\nu {\rm (erg/s/Hz)}$"
        y, yerr = compute_lnu_error(intensity,errors)
    elif (yunit == 'counts'):
        ylabel = r"$N_\nu {\rm (counts/s/Hz)}$"
        y, yerr = compute_counts_error(intensity,nu,errors)
    elif (yunit == 'polfrac'):
        ylabel = r"$\rm Pol.\; Fraction \; (\%)$"
        y, yerr = compute_pol_frac_error(intensity,errors)
    elif (yunit == 'polangle'):
        ylabel = r"$\rm Pol.\; Angle$"
        y, yerr = compute_pol_angle_error(intensity,errors)
    elif (yunit == 'q'):
        ylabel = r"$Q_\nu/I_\nu$"
        y, yerr = compute_q_error(intensity,errors)
    elif (yunit == 'u'):
        ylabel = r"$U_\nu/I_\nu$"
        y, yerr = compute_u_error(intensity,errors)
    elif (yunit == 'fluxfrac'):
        ylabel = r"$I_\nu/F_\nu$"
        y, yerr = compute_flux_frac_error(intensity,xfaces,errors)

    # Return x and y variables, there labels, and possible error on y
    return x,y,yerr,xlabel,ylabel


def make_plot(x,y,yerr=None,ax=None,xmin=None,xmax=None,ymin=None,ymax=None,
              xlabel=None,ylabel=None,xscale=None,yscale=None,fmt=None,
              **kwargs):
    """
    General wrapper for plotting of Monte Carlo spectral plots
    """

    if (ax is None):
        # Create figure, axis and assume a single plot window
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if (fmt is None):
        fmt = '.'

    if (yerr is not None):
        ax.errorbar(x,y,yerr=yerr,fmt=fmt,**kwargs)
    else:
        ax.plot(x,y,fmt,**kwargs)

    # Set axis labelx
    if (xlabel is not None):
        ax.set_xlabel(xlabel)
    if (ylabel is not None):
        ax.set_ylabel(ylabel)

    # Set axis scales
    if (xscale is not None):
        ax.set_xscale(xscale)
    if (yscale is not None):
        ax.set_yscale(yscale)

    # (re)Set plot ranges
    if xmin is not None:
        if xmax is not None:
            ax.set_xlim(xmin,xmax)
        else:
            ax.set_xlim(left=xmin)
    elif xmax is not None:
        ax.set_xlim(right=xmax)

    if ymin is not None:
        if ymax is not None:
            ax.set_ylim(ymin,ymax)
        else:
            print(ymin)
            ax.set_ylim(bottom=ymin)
    elif ymax is not None:
        ax.set_ylim(top=ymax)

    return ax


def get_luminosity(spec):
    """
    Computes the integrated luminosity corresponding to a spectrum
    """

    nx = spec['nx']
    nmu = spec['nmu']
    nphi = spec['nphi']
    mumid = 0.5*(spec['mufaces'][1:]+spec['mufaces'][:-1])

    # Compute frequency width and mean energy (in ergs) of bins
    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    xaxis = spec['units']
    xfaces = spec['xfaces']
    if (xaxis == 'kev'):
        dnu = (xfaces[1:]-xfaces[:-1])*1000.*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*1000.*everg
    if (xaxis == 'ev'):
        dnu = (xfaces[1:]-xfaces[:-1])*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*everg
    if (xaxis == 'nu'):
        dnu = (xfaces[1:]-xfaces[:-1])
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*h
    if (xaxis == 'lambda'):
        dnu = (1./xfaces[:-1]-1./xfaces[1:])*c/1.e8
        emid = 0.5*(1./xfaces[:-1]-1./xfaces[1:])*c*h/1.e8

    # compute sum over frequency and solid angle
    lumin = 0.
    for k in range(nx):
        for j in range(nmu):
            for i in range(nphi):
                lumin += dnu[k]*mumid[j]*spec['intensity'][0,i,j,k]
    lumin *= 2.*np.pi/nmu/nphi
    return lumin

def build_bins(xmin,xmax,nx,logx):
    """
    Builds a x-axis grid for binning the photons
    """
    if (logx):
        return np.logspace(np.log10(xmin),np.log10(xmax),nx+1)
    else:
        return np.linspace(xmin,xmax,nx+1)

def get_bins(xphots,xfaces,nx,uniform=True,log=True):

    xbins = np.zeros(len(xphots),dtype=int)-1
    if uniform:
        if log:
            xlfaces = np.log10(xfaces)
            xwidth = xlfaces[nx]-xlfaces[0]
            for i,xphot in enumerate(xphots):
                xbins[i] = int((np.log10(xphot)-xlfaces[0])/xwidth*float(nx))
                if ((xbins[i] < 0) or (xbins[i] >= nx)):
                    xbins[i] = -1
        else:
            xwidth = xfaces[nx]-xfaces[0]
            for i,xphot in enumerate(xphots):
                if np.isinf(xphot):
                    xbins[i] = -1
                else:
                    xbins[i] = int((xphot-xfaces[0])/xwidth*float(nx))
                if ((xbins[i] < 0) or (xbins[i] >= nx)):
                    xbins[i] = -1
        return xbins
    else:
        # use binary search for non uniform data
        return get_bins_binary_search(xphots,xfaces,nx)

def get_bins_binary_search(xphots,xfaces,nx):
    """
    Returns x bin numbers corresponding to xphots for non-uniformly
    binned data via binary search.
    """
    # Exclude values outside of search range
    indsp = (xphots > xfaces[nx]).nonzero()
    indsm = (xphots < xfaces[0]).nonzero()
    xbins = np.searchsorted(xfaces,xphots)-1
    xbins[indsp] = -1
    xbins[indsm] = -1
    return xbins

def get_angle_bins_cartesian(photons,nmu,mufaces,nphi,phifaces):
    """
    Bin angles in theta, phi defined relative to x,y,z axes
    """

    if ((nmu == 1) and (mufaces[0] <= 0.) and (mufaces[1] >= 1.0)):
        skipmu = True
    else:
        skipmu = False
    if ((nphi == 1) and (phifaces[0] <= 0.) and (phifaces[1] >= 2.*np.pi)):
        skipphi = True
    else:
        skipphi = False

    if (skipmu and skipphi):
        return np.zeros(photons.nphot,dtype=int),np.zeros(photons.nphot,dtype=int)

    if photons.coord == 'spherical_polar':
        kr = photons.k1
        kth = photons.k2
        kph = photons.k3
        cth = np.cos(photons.x2)
        sth = np.sin(photons.x2)
        cph = np.cos(photons.x3)
        sph = np.sin(photons.x3)
        if (not skipphi):
            kx = kr*sth*cph + kth*cth*cph - kph*sph
            ky = kr*sth*sph + kth*cth*sph - kph*cph
        if (not skipmu):
            kz = kr*cth - kth*sth
    else:
        if (not skipphi):
            kx = photons.k1
            ky = photons.k2
        if (not skipmu):
            kz = photons.k3

    if (skipmu):
        # return 0
        mubins = np.zeros(photons.nphot,dtype=int)
    else:
        # Bin based on k . z
        mu = abs(kz)
        mubins = get_bins(mu,mufaces,nmu,log=False)
    if (skipphi):
        # return 0
        phibins = np.zeros(photons.nphot,dtype=int)
    else:
        if (skipmu):
            mu = abs(kz)
        smu = np.sqrt(1.-mu*mu)

        phi = np.arctan2(ky,kx)
        phi[(phi<0.).nonzero()] += 2.*np.pi
        phibins = get_bins(phi,phifaces,nphi,log=False)


    return mubins,phibins


def get_angle_bins_spherical(photons,nmu,mufaces):
    """
    Bin angles relative to local radial direction.  Here we only bin
    in polar angle.
    """
    if ((nmu == 1) and (mufaces[0] <= 0.) and (mufaces[1] >= 1.0)):
        skipmu = True
    else:
        skipmu = False

    if (skipmu):
        return np.zeros(photons.nphot,dtype=int),np.zeros(photons.nphot,dtype=int)

    if photons.coord == 'spherical_polar':
        kr = photons.k1
    else:
        cth = np.cos(photons.x2)
        sth = np.sin(photons.x2)
        cph = np.cos(photons.x3)
        sph = np.sin(photons.x3)
        kr = sth*(cph*photons.k1+sph*photons.k2)+cth*photons.k3

    # Bin based on k_r
    mu = kr
    mubins = get_bins(mu,mufaces,nmu,log=False)

    # return 0 for phi
    phibins = np.zeros(photons.nphot,dtype=int)

    return mubins,phibins

def get_angle_bins_hybrid(photons,nmu,mufaces,nphi,phifaces):
    """
    Bin angles in theta, phi defined relative to x,y,z axes, but with phi
    defined by local azimuthal angle
    """
    if (nmu == 1):
        printf("Error: this function requires nmu > 1")
        return  np.zeros(photons.nphot,dtype=int),np.zeros(photons.nphot,dtype=int)
    if (nphi == 1):
        printf("Error: this function requires nphi > 1")
        return np.zeros(photons.nphot,dtype=int),np.zeros(photons.nphot,dtype=int)

    if (photons.coord != 'spherical_polar'):
        printf("Error: this function only works with spherical polar")
        return np.zeros(photons.nphot,dtype=int),np.zeros(photons.nphot,dtype=int)
    else:
        kr = photons.k1
        kth = photons.k2
        kph = photons.k3
        cth = np.cos(photons.x2)
        sth = np.sin(photons.x2)
        kz = kr*cth - kth*sth

    # Bin based on k . z
    mu = abs(kz)
    mubins = get_bins(mu,mufaces,nmu,log=False)

    phi = np.arcsin(kph)
    # inds should be empty for outgoing photons
    inds = (kr < 0.).nonzero()
    phi[inds] = np.pi - phi[inds]
    phi[(phi<0.).nonzero()] += 2.*np.pi
    phibins = get_bins(phi,phifaces,nphi,log=False)

    return mubins,phibins

def make_spectrum(phots,nx,xmin,xmax,xaxis='kev',logx=True,nmu=1,mumin=0,mumax=1.,
                  nphi=1,phimin=0,phimax=2.*np.pi,yerror=True,mask=None,
                  xfunc=None,anglebin='cartesian',**kwargs):
    """
    Makes spectrum (dict) from photon object
    """

    # Store spectrum as a dictionary
    spectrum = {}

    # Store integration time
    spectrum['dt'] = phots.dt

    # Store total number of photons for refernce
    spectrum['ntot'] = phots.ntot

    # Set x binning variable
    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    preset = False
    spectrum['xaxis'] = xaxis
    if (xaxis == 'kev'):
        xphots = phots.energy/everg/1000.
        preset = True
    if (xaxis == 'ev'):
        xphots = phots.energy/everg
        preset = True
    if (xaxis == 'nu'):
        xphots = phots.energy/h
        preset = True
    if (xaxis == 'lambda'):
        xphots = c*h/(phots.energy*1.e8)
        preset = True
    if (not preset):
        if xfunc is None:
            raise RuntimeError('Unrecognized xunit and no xfunc provided')
        else:
            xphots = xfunc(phots.energy,True,**kwargs)

    # Create bins
    xfaces = build_bins(xmin,xmax,nx,logx)
    spectrum['nx'] = nx
    spectrum['xfaces'] = xfaces

    # Get x bins
    xbins = get_bins(xphots,xfaces,nx,log=logx)

    # Get angle bins
    spectrum['nmu'] = nmu
    mufaces = build_bins(mumin,mumax,nmu,False)
    spectrum['mufaces'] = mufaces

    spectrum['nphi'] = nphi
    phifaces = build_bins(phimin,phimax,nphi,False)
    spectrum['phifaces'] = phifaces

    if (anglebin == 'cartesian'):
        mubins, phibins = get_angle_bins_cartesian(phots,nmu,mufaces,nphi,phifaces)
    elif (anglebin == 'spherical'):
        mubins, phibins = get_angle_bins_spherical(phots,nmu,mufaces)
    elif (anglebin == 'hybrid'):
        mubins, phibins = get_angle_bins_hybrid(phots,nmu,mufaces,nphi,phifaces)
    else:
        print("Error: anglebin == "+anglebin+". Must be cartesian, spherical, or hybrid.")

    # Create intensity grid and loop over photons to add contribution
    nintens = 1
    if phots.polarized:
        spectrum['polarized'] = 'true'
        nintens += 2
    else:
        spectrum['polarized'] = 'false'

    spectrum['nintens'] = nintens
    count = np.zeros((nphi,nmu,nx))
    intensity = np.zeros((nintens,nphi,nmu,nx))
    if yerror:
        errors = np.zeros((nintens,nphi,nmu,nx))

    if (mask is not None):
        xbins[mask] = -1

    for i in range(phots.nphot):
        if ((xbins[i] >= 0) and (mubins[i] >= 0) and (phibins[i] >= 0)):
            wght = phots.weight[i]
            if ((phots.q[i]**2+phots.u[i]**2) > 1.001):
                wght = 0.
                print("Warning: polarization too high: ",phots.q[i],phots.u[i],phots.weight[i],np.sqrt(phots.q[i]**2+phots.u[i]**2))
            #print phibins[i],mubins[i],xbins[i]
            count[phibins[i],mubins[i],xbins[i]] += 1.
            intensity[0,phibins[i],mubins[i],xbins[i]] += wght
            if phots.polarized:
                intensity[1,phibins[i],mubins[i],xbins[i]] += wght*phots.q[i]
                intensity[2,phibins[i],mubins[i],xbins[i]] += wght*phots.u[i]
            if yerror:
                errors[0,phibins[i],mubins[i],xbins[i]] += (wght)**2
                if phots.polarized:
                    errors[1,phibins[i],mubins[i],xbins[i]] += (wght*phots.q[i])**2
                    errors[2,phibins[i],mubins[i],xbins[i]] += (wght*phots.u[i])**2

    # Compute frequency width and mean energy (in erg) of bins
    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    if (xaxis == 'kev'):
        dnu = (xfaces[1:]-xfaces[:-1])*1000.*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*1000.*everg
    if (xaxis == 'ev'):
        dnu = (xfaces[1:]-xfaces[:-1])*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*everg
    if (xaxis == 'nu'):
        dnu = (xfaces[1:]-xfaces[:-1])
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*h
    if (xaxis == 'lambda'):
        dnu = (1./xfaces[:-1]-1./xfaces[1:])*c/1.e8
        emid = 0.5*(1./xfaces[:-1]-1./xfaces[1:])*c*h/1.e8
    if (not preset):
        efaces = xfunc(xfaces,False,**kwargs)
        emid = 0.5*(efaces[1:]+efaces[:-1])
        dnu  = (efaces[1:]-efaces[:-1])/h

    # Normalize intensities
    mumid = 0.5*(mufaces[1:]+mufaces[:-1])
    for k in range(nintens):
        for j in range(nphi):
            for i in range(nmu):
                fac = nphi*nmu*emid/(mumid[i]*dnu*2.*np.pi*phots.dt)
                intensity[k,j,i,:] *= fac
                if yerror:
                    errors[k,j,i,:] *= fac**2
    spectrum['intensity'] = intensity

    # Finish computing errors on intensities
    if yerror:
        for k in range(nintens):
            inds = np.where(count > 1.)
            errors[k,inds[0],inds[1],inds[2]] = np.sqrt(errors[k,inds[0],inds[1],inds[2]]
              - intensity[k,inds[0],inds[1],inds[2]]**2 / count[inds[0],inds[1],inds[2]])
            inds = np.where(count <= 1.)
            errors[k,inds[0],inds[1],inds[2]] = 0.
        spectrum['errors'] = errors

    if yerror:
        spectrum['yerror'] = "true"
    else:
        spectrum['yerror'] = "false"

    return spectrum

def get_image_bins(phots,rcam,ifaces,xfaces,yfaces):
    """
    Bin photons in image plane coordinates
    """
    ninc = ifaces.size - 1
    nx  = xfaces.size - 1
    ny = yfaces.size - 1

    thc = 0.5*(xfaces[1:]+xfaces[:-1])
    if (phots.coord == 'spherical_polar'):
        sth = np.sin(phots.x2)
        cth = np.cos(phots.x2)
        sph = np.sin(phots.x3)
        cph = np.cos(phots.x3)
        xp = phots.x1*sth*cph
        yp = phots.x1*sth*sph
        zp = phots.x1*cth
        rp = phots.x1
        kdx = rp*phots.k1
        kx = phots.k1*sth*cph + phots.k2*cth*cph - phots.k3*sph
        ky = phots.k1*sth*sph + phots.k2*cth*sph + phots.k3*cph
        kz = phots.k1*cth - phots.k2*sth
    elif (phots.coord == 'cartesian'):
        xp = phots.x1
        yp = phots.x2
        zp = phots.x3
        rp = np.sqrt(xp**2+yp**2+zp**2)
        kx = phots.k1
        ky = phots.k2
        kz = phots.k3
        kdx = xp*kx+yp*ky+zp*kz

    dl = np.sqrt(rcam**2-rp**2+kdx**2)-kdx
    xf = xp + dl * kx
    yf = yp + dl * ky
    zf = zp + dl * kz

    cthf = zf/rcam
    sthf = np.sqrt(1.-cthf**2)
    phf =  np.arctan2(yf,xf)
    cphf = np.cos(phf)
    sphf = np.sin(phf)
    #np.set_printoptions(threshold=1000)

    ibins = get_bins(cthf,ifaces,ninc,log=False)

    kth = kx*cthf*cphf + ky*cthf*sphf - kz*sthf
    kph = -kx*sphf + ky*cphf
    norm = np.sqrt(1.+kth**2+kph**2)
    y = kth*rcam*norm
    x = kph*rcam*norm

    xbins = get_bins(x,xfaces,nx,log=False)
    ybins = get_bins(y,yfaces,ny,log=False)
    #for i,q in enumerate(cthf):
    #    print(q,ibins[i],x[i],y[i],xbins[i],ybins[i])
    return ibins, xbins, ybins

def make_image_mc(phots,rcam,ninc,imin,imax,nen,emin,emax,
                  nx,xmin,xmax,ny,ymin,ymax,unit,mask=None,**kwargs):
    """
    Create a binned image from photon list
    """

    # Store the image as a dictionary
    image = {}

    # Store integration time
    image['dt'] = phots.dt

    # Store total number of photons for refernce
    image['ntot'] = phots.ntot

    # Create bins for viewer inclination
    ifaces = build_bins(imin,imax,ninc,False)
    image['ninc'] = ninc
    image['ifaces'] = ifaces

    # Create bins for observed frequency/photon energy
    efaces = build_bins(emin,emax,nen,True)
    image['nen'] = nen
    image['efaces'] = efaces

    # Create bins for image plane, image will be uniform 2d array
    image['nx'] = nx
    image['ny'] = ny
    xfaces = build_bins(xmin,xmax,nx,False)
    yfaces = build_bins(ymin,ymax,ny,False)
    image['xfaces'] = xfaces
    image['yfaces'] = yfaces

    # set units
    image['unit'] = unit

    ibins, xbins, ybins = get_image_bins(phots,rcam,ifaces,xfaces,yfaces)
    #set ebins temporarily to 0
    ebins = np.zeros(len(ibins),dtype=int)

    # Create intensity grid and loop over photons to add contribution
    nintens = 1
    if phots.polarized:
        image['polarized'] = True
        nintens += 2
    else:
        image['polarized'] = False
    image['nintens'] = nintens


    if (mask is not None):
        ibins[mask] = -1

    intensity = np.zeros((nintens,ninc,nen,ny,nx))
    for i in range(phots.nphot):
        if ((ibins[i] >= 0) and (ebins[i] >= 0) and (xbins[i] >= 0) and (ybins[i] >= 0)):
            # Weight includes energy -- slightly different from spectra
            wght = phots.weight[i]*phots.energy[i]
            intensity[0,ibins[i],ebins[i],ybins[i],xbins[i]] += wght
            if phots.polarized:
                intensity[1,ibins[i],ebins[i],ybins[i],xbins[i]] += wght*phots.q[i]
                intensity[2,ibins[i],ebins[i],ybins[i],xbins[i]] += wght*phots.u[i]

    # Normalize intensities
    mumid = abs(0.5*(ifaces[1:]+ifaces[:-1]))
    dmu = ifaces[1:]-ifaces[:-1]

    if (image['nen'] == 1):
        dnu = np.array([1.])
    else:
        h = 6.62607015e-27
        everg = 1.6021772e-12
        dnu = (efaces[1:]-efaces[:-1])*1000.*everg/h

    dx = xfaces[1:]-xfaces[:-1]
    dy = yfaces[1:]-yfaces[:-1]
    area = np.outer(dy,dx)
    for k in range(nintens):
        for j in range(ninc):
            for i in range(nen):
                # not divided by dnu for now
                fac = dnu[i]*dmu[j]*mumid[j]*2.*np.pi*phots.dt
                intensity[k,j,i,:,:] /= fac*area


    image['intensity'] = intensity

    return image

def subsample_polarization(q,u,x,y,step,average):
    """
    Subsampling of itensity array
    """

    nx = len(x)
    ny = len(y)

    if (not average):
        x = x[step // 2:nx:step]
        y = y[step // 2:ny:step]
        q = q[step // 2:nx:step,step // 2:ny:step]
        u = u[step // 2:nx:step,step // 2:ny:step]
        return q,u,x,y
    else:
        # too lazy to work out pythony way of doing this
        xp = np.zeros(nx // step)
        yp = np.zeros(ny // step)
        qp = np.zeros((nx // step,ny // step))
        up = np.zeros((nx // step,ny // step))
        for i in range(nx // step):
            xp[i] = np.average(x[i*step:(i+1)*step])
        for i in range(ny // step):
            yp[i] = np.average(y[i*step:(i+1)*step])
        for i in range(nx // step):
            for j in range(ny // step):
                qp[i,j] = np.average(q[i*step:(i+1)*step,j*step:(j+1)*step])
                up[i,j] = np.average(u[i*step:(i+1)*step,j*step:(j+1)*step])

        return qp,up,xp,yp



def plot_image(image,iinc,type='intensity',pvec=False,average=False,step=4,
               ax=None,**kwargs):
    """
    Plot an image
    """
    if (ax is None):
        # Create figure, axis and assume a single plot window
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    vmin = kwargs['vmin']
    vmax = kwargs['vmax']
    cmap = plt.get_cmap(kwargs['colormap'])
    plt.figure()

    if (polarization_requested(type)):
        if (not image['polarized']):
            raise RuntimeError("Polarization type requested ("+type+
                               ") but image is unpolarized")
    if (type == 'intensity'):
        vals = image['intensity'][0,iinc,0,:,:]
        clabel=r"$I$"
        if vmin is None:
            vmin = 1.e-5*np.max(vals)
    elif (type == 'q'):
        vals = image['intensity'][1,iinc,0,:,:]
        clabel=r"$Q/I$"
    elif (type == 'u'):
        vals = image['intensity'][2,iinc,0,:,:]
        clabel=r"$U/I$"
    elif (type == 'polangle'):
        q = image['intensity'][1,iinc,0,:,:]
        u = image['intensity'][2,iinc,0,:,:]
        vals = 90./np.pi*np.arctan2(u,q)
        vals[vals < 0.] += 360.
        if (vmin is None):
            vmin = 0.
        clabel = r"$\rm Pol.\; Angle$"
    elif (type == 'polfrac'):
        q = image['intensity'][1,iinc,0,:,:]
        u = image['intensity'][2,iinc,0,:,:]
        vals = np.sqrt(q**2+u**2)
        if (vmin is None):
            vmin = 0.
        clabel = r"$\rm Pol.\; Fraction$"
    else:
        raise RuntimeError("Type:"+type+" is not defined.")

    if kwargs['vnorm']:
        vm = np.max(vals)
        vals = vals / vm
        vmax = 1.
        vmin /= vm

    if (kwargs['logc']):
        vals[vals <= 0.] = 1.e-20 * np.max(vals)
        norm = colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin,vmax=vmax)

    x = 0.5*(image['xfaces'][1:]+image['xfaces'][:-1])
    y = 0.5*(image['yfaces'][1:]+image['yfaces'][:-1])
    x_2d, y_2d = np.meshgrid(x,y)
    im = plt.pcolormesh(x_2d, y_2d, vals, cmap=cmap, norm=norm)

    plt.xlim(image['xfaces'][0],image['xfaces'][-1])
    plt.ylim(image['yfaces'][0],image['yfaces'][-1])
    if (image['unit'] == 'cm'):
        plt.xlabel(r"$x \; (\rm cm)$")
        plt.ylabel(r"$y \; (\rm cm)$")
        if (type == 'intensity'):
            clabel=r"$I \; (\rm erg/s/cm^2)$"
    else:
        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")

    plt.colorbar(im,label=clabel)
    plt.gca().set_aspect('equal')
    if (pvec):
        if (image['polarized']):
            q = image['intensity'][1,iinc,0,:,:]
            u = image['intensity'][2,iinc,0,:,:]
            q, u, x, y = subsample_polarization(q,u,x,y,step,average)
            x_pol, y_pol = np.meshgrid(x,y)

            pol_angle = 0.5 * np.arctan2(u,q)
            pol_frac = np.sqrt(q*q+u*u)
            vx = pol_frac*np.cos(pol_angle)
            vy = pol_frac*np.sin(pol_angle)

            plt.quiver(x_pol, y_pol, vx, vy, color='k',headwidth=0, headlength=0,
                       headaxislength=0, scale = None,pivot='middle')
        else:
            raise RuntimeError("Polarization vectors requested but image is unpolarized")

def write_image(filename,image):
    """
    Writes image to output file
    """

    # Open outfile
    outfile = open(filename, 'w')

    ninc = image['ninc']
    nen = image['nen']
    nx = image['nx']
    ny = image['ny']

    # Write header information
    outfile.write("dt={:.8e}\n".format(image['dt']))
    outfile.write("ninc={:d}\n".format(ninc))
    outfile.write("nen={:d}\n".format(nen))
    outfile.write("nx={:d}\n".format(nx))
    outfile.write("ny={:d}\n".format(ny))
    outfile.write("unit="+image['unit']+"\n")
    outfile.write("ntot={:d}\n".format(image['ntot']))
    outfile.write("nintens={:d}\n".format(image['nintens']))
    if image['polarized']:
        outfile.write("polarized=true\n")
    else:
        outfile.write("polarized=false\n")
    outfile.close()

    # Write binfaces
    outfile = open(filename, 'ab')
    myfmt='>'+'d'*(ninc+1)
    bin=struct.pack(myfmt,*(image['ifaces']))
    outfile.write(bin)
    myfmt='>'+'d'*(nen+1)
    bin=struct.pack(myfmt,*(image['efaces']))
    outfile.write(bin)
    myfmt='>'+'d'*(nx+1)
    bin=struct.pack(myfmt,*(image['xfaces']))
    outfile.write(bin)
    myfmt='>'+'d'*(ny+1)
    bin=struct.pack(myfmt,*(image['yfaces']))
    outfile.write(bin)
    # Write data
    nelements = (image['nintens']*ninc*nen*ny*nx)
    myfmt='>'+'d'*nelements
    bin=struct.pack(myfmt,*(image['intensity'].reshape(nelements)))
    outfile.write(bin)
    outfile.close()

def read_image(filename):
    """
    Read image and return as a dictionary
    """

    # Read raw data
    with open(filename, 'rb') as data_file:
        raw_data = data_file.read()
    raw_data_ascii = raw_data.decode('ascii', 'replace')

    image = {}
    current_index = 0

    # Function for skipping though the file
    def skip_string(expected_string):
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index+expected_string_len] != expected_string:
            raise RuntimeError('File not formatted as expected')
        return current_index+expected_string_len

    try:
        current_index = skip_string("dt=")
        end_of_line_index = current_index + 1
        while raw_data_ascii[end_of_line_index] != '\n':
            end_of_line_index += 1
        image['dt'] = list(map(float,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
        current_index = end_of_line_index + 1
    except:
        print("Image file contains no dt entry. Setting to 1.")
        image['dt'] = 1.

    current_index = skip_string("ninc=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['ninc'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1

    current_index = skip_string("nen=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['nen'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1

    current_index = skip_string("nx=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['nx'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1

    current_index = skip_string("ny=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['ny'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1

    current_index = skip_string("unit=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['unit'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1

    current_index = skip_string("ntot=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['ntot'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1

    current_index = skip_string("nintens=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['nintens'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1

    current_index = skip_string("polarized=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    image['polarized'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1
    if image['polarized'] == 'false':
        image['polarized'] = False
    if image['polarized'] == 'true':
        image['polarized'] = True

    # Read in faces
    ninc = image['ninc']
    format_string = '>' + 'd'*(ninc+1)
    begin_index = current_index
    end_index = begin_index + 8*(ninc+1)
    image['ifaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    nen = image['nen']
    format_string = '>' + 'd'*(nen+1)
    begin_index = end_index
    end_index = begin_index + 8*(nen+1)
    image['efaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    nx = image['nx']
    format_string = '>' + 'd'*(nx+1)
    begin_index = end_index
    end_index = begin_index + 8*(nx+1)
    image['xfaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    ny = image['ny']
    format_string = '>' + 'd'*(ny+1)
    begin_index = end_index
    end_index = begin_index + 8*(ny+1)
    image['yfaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))

    # Read intensities
    nintens = image['nintens']
    nelements = nintens*ninc*nen*ny*nx
    format_string = '>' + 'd'*nelements
    begin_index = end_index
    end_index = begin_index + 8*nelements
    vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    image['intensity'] = vals.reshape((nintens,ninc,nen,ny,nx))
    return image
