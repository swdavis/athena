"""
Support for manipulating Monce Carlo photon trajectory lists
"""

# standard python modules
import numpy as np
import struct
import matplotlib.pyplot as plt

def read_list(filename):
    """
    Read unformated list output and return as a dictionary
    """
    # Read raw data
    with open(filename, 'rb') as data_file:
        raw_data = data_file.read()
    raw_data_ascii = raw_data.decode('ascii', 'replace')

    # Store in dictionary
    phtraj = {}
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
    phtraj['length'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("maxstep=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phtraj['maxstep'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("npars=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phtraj['npars'] = map(int,raw_data_ascii[current_index:end_of_line_index].split(' '))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("coord=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    phtraj['coord'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1
    # Read in steps
    format_string = '>' + 'i'*phtraj['length']
    begin_index = current_index
    end_index = begin_index + 4 *phtraj['length']
    nsteps = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]),dtype=int)
    phtraj['nsteps'] = nsteps
    current_index = end_index
    nsteptot = 0
    for i in range(phtraj['length']):
        nsteptot += nsteps[i]
    # Read in data
    nelements = nsteptot * phtraj['npars']
    format_string = '>' + 'd'*nelements
    begin_index = current_index
    end_index = begin_index + 8*nelements
    vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    phtraj['list'] = np.zeros((phtraj['length'],phtraj['maxstep'],phtraj['npars']))
    n = 0
    for i in range(phtraj['length']):
        for j in range(nsteps[i]):
            for k in range(phtraj['npars']):
                phtraj['list'][i,j,k] = vals[n]
                n += 1

    return phtraj

def write_list(filename,phtraj):
    """
    Write photon list (dictionary) to file
    """
    outfile = open(filename, 'w')
    
    # Write header information
    outfile.write("length={:d}\n".format(phtraj['length']))
    outfile.write("maxstep={:d}\n".format(phtraj['maxstep']))
    outfile.write("npars={:d}\n".format(phtraj['npars']))
    outfile.write("coord="+phtraj['coord']+"\n")
    # Compute total steps
    nsteptot = 0
    for i in range(phtraj['length']):
        nsteptot += phtraj['nsteps'][i]
    nelements = nsteptot*phtraj['npars']
    # Write number of steps
    myfmt = '>' + 'i'*phtraj['length']
    bin=struct.pack(myfmt,*(phtraj['nsteps']))
    outfile.write(bin)
    # Write list data
    traj = phtraj['list'] = phtraj['list'].reshape(nelements)
    myfmt='>'+'d'*nelements
    bin=struct.pack(myfmt,*(traj))
    outfile.write(bin)
    outfile.close()

def plot_trajectory_projection(traj,axs=None,l0=1.):

    if (axs is None):
        # Create figure, axes
        fig, axs = plt.subplots(1, 3, figsize=(14,4))
    print traj['list'].shape
    for i in range(traj['length']):
        n = traj['nsteps'][i]
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        if (traj['coord'] == 'cartesian'):
            for j in range(n):
                x[j] = traj['list'][i,j,0]
                y[j] = traj['list'][i,j,1]
                z[j] = traj['list'][i,j,2]
        if (traj['coord'] == 'spherical_polar'):
            for j in range(n):
                r = traj['list'][i,j,0]
                cth = np.cos(traj['list'][i,j,1])
                sth = np.sin(traj['list'][i,j,1])            
                cph = np.cos(traj['list'][i,j,2])
                sph = np.sin(traj['list'][i,j,2])
                x[j] = r*sth*cph
                y[j] = r*sth*sph
                z[j] = r*cth
        if (traj['coord'] == 'kerr-schild'):
            for j in range(n):
                r = traj['list'][i,j,0]
                cth = np.cos(traj['list'][i,j,1])
                sth = np.sin(traj['list'][i,j,1])            
                cph = np.cos(traj['list'][i,j,2])
                sph = np.sin(traj['list'][i,j,2])
                x[j] = r*sth*cph
                y[j] = r*sth*sph
                z[j] = r*cth

        axs[0].plot(x,y,'r.', markersize=1.0)
        axs[1].plot(x,z,'r.', markersize=1.0)
        axs[2].plot(y,z,'r.', markersize=1.0)
        
        axs[0].set_xlim([-l0,l0])
        axs[1].set_xlim([-l0,l0])
        axs[2].set_xlim([-l0,l0])
        axs[0].set_ylim([-l0,l0])
        axs[1].set_ylim([-l0,l0])
        axs[2].set_ylim([-l0,l0])

        #plt.show()
        return axs
