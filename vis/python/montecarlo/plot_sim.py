#! /usr/bin/env python

"""
Plot athena++ output variables in multiple panels
"""

# python standard modules
import sys
import argparse
import numpy as np
import string as str
import math
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.colors as colors

# athena++ specific
import athena_read

def par_get_string(par,block,pardef,input):
   """
   Parse input for string parameter. If found set parameter, else set default.
   """
   # Simply a wrapper for par_get()
   val = par_get(par,block,input)
   if ((val == 'ERROR') or (val == 'NOTFOUND')):
      if (pardef == 'NODEF'):            
         print("Required parameter '{:s}' not found in input.".format(par))
         sys.exit(0)
      else:
         return pardef
   else:
      return val


def par_get_float(par,block,pardef,input):
   """
   Parse input for float parameter. If found set parameter, else set default.
   """
   val = par_get(par,block,input)
   if ((val == 'ERROR') or (val == 'NOTFOUND')):
      if (pardef == 'NODEF'):            
         print("Required parameter '{:s}' not found in input.".format(par))
         sys.exit(0)
      else:
         return pardef
   else:
      try:
         return float(val)
      except:
         print("Invalid float: '{:s}'. Returning default.".format(val))
         return pardef


def par_get_int(par,block,pardef,input):
   """
   Parse input for integer parameter. If found set parameter, else set default.
   """
   val = par_get(par,block,input)
   if ((val == 'ERROR') or (val == 'NOTFOUND')):
      if (pardef == 'NODEF'):            
         print("Required parameter '{:s}' not found in input.".format(par))
         sys.exit(0)
      else:
         return pardef
   else:
      try:
         return int(val)
      except:
         print("Invalid integer: '{:s}'. Returning default.".format(val))
         return pardef
         
def par_get_boolean(par,block,pardef,input):
   """
   Parse input for boolean parameter. If found set parameter, else set default.
   """
   val = par_get(par,block,input)
   if ((val == 'ERROR') or (val == 'NOTFOUND')):
      if (pardef == 'NODEF'):            
         print("Required parameter '{:s}' not found in input.".format(par))
         sys.exit(0)
      else:
         return pardef
   else:
      if (val == "True"):
         return True
      elif (val == "true"):
         return True
      elif (val == "1"):
         return True
      elif (val == "False"):
         return False
      elif (val == "false"):
         return false
      elif (val == "0"):
         return false
      else:
         print("Parameter '{:s}' assigned non-boolean value: " \
               "{:s}. Returning default".format(par,val) )
         return pardef

def par_get(par,block,input):
   """
   Parse input for parameter. If found return parameter, else return flag.
   """

   # First determine if in block
   found_par = False
   in_block = False
   for line in input:
      if line[0] == '#':
         # ignore comments
         continue
      if "<"+block+">" in line:
         # start of target block
         in_block = True
         continue
      if in_block:
         if "<" in line:
            # start of different block
            in_block =False
            continue
         else:
            if par in line:
               #if (par == str.split(line)[0]):
               if (par == line.split()[0]):
                  found_par = True
                  break

   if (found_par):
      list = line.split()
      if (len(list) < 3):
         print("Syntax error: '{:2}' Returning default.".format(line.strip()))
         return 'ERROR'
      elif (list[1] != "=" ):  
         print("Syntax error: '{:2}' Returning default.".format(line.strip()))
         return 'ERROR'
      else:
         return list[2]
   else:
      return 'NOTFOUND'


def set_plot_geometry(infname, args):
   """
   Parses input file if given and sets plot geometery
   """

   # create empty dictioary
   argsn = {}

   if (infname is None):
      input = ""
   else:
      # Read the input file into an array of strings -- one line for each element
      infile = open(infname,'r')
      input = infile.readlines()
      infile.close()

   # Determine numbers of columns and rows
   argsn['npanels'] = npanels = par_get_int('npanels','geometry',args.pop('npanels'),input)
   argsn['ncol'] = ncol = par_get_int('ncol','geometry',args.pop('ncol'),input)
   if (ncol is None):
      ncol = math.isqrt(npanels)
   nrow = npanels // ncol
   if (npanels % ncol > 0):
      nrow += 1
   argsn['nrow'] = nrow

   # Set figure height and width
   argsn['figheight'] = par_get_float('figheight','geometry',args.pop('figheight'),input)
   argsn['figwidth'] = par_get_float('figwidth','geometry',args.pop('figwidth'),input)
   if (argsn['figwidth'] is None):
      if (argsn['figheight'] is None):
         argsn['figwidth'] = min(ncol * 4.,12.)
         argsn['figheight'] = min(nrow * 3.8,11.4)
      else:
         argsn['figwidth'] = min(argsn['figheight'] * nrow / ncol * 4.,12.)

   aspect = min(argsn['figwidth']/nrow,argsn['figwidth']/ncol)
   deftextsize = int(2.*aspect)
   
   argsn['textsize'] = par_get_int('textsize','geometry',args.pop('textsize'),input)
   if argsn['textsize'] is None:
      argsn['textsize'] = deftextsize

   argsn['axes_all'] = par_get_float('axes_all','geometry',args.pop('axes_all'),input)
   argsn['hspace'] = par_get_float('hspace','geometry',args.pop('hspace'),input)
   argsn['wspace'] = par_get_float('wspace','geometry',args.pop('wspace'),input)
   argsn['top'] = par_get_float('top','geometry',args.pop('top'),input)
   argsn['bottom'] = par_get_float('bottom','geometry',args.pop('bottom'),input)
   argsn['right'] = par_get_float('right','geometry',args.pop('right'),input)
   argsn['left'] = par_get_float('left','geometry',args.pop('left'),input)

   if (argsn['hspace'] is None):
      if (argsn['axes_all']):
         argsn['hspace'] = 0.15
      else:
         argsn['hspace'] = 0.05
   if (argsn['wspace'] is None):
      if (argsn['axes_all']):
         argsn['wspace'] = 0.5
      else:
         argsn['wspace'] = 0.3

   if (argsn['top'] is None):
      argsn['top'] = 0.95
   if (argsn['bottom'] is None):
      if (nrow == 1):
         argsn['bottom'] = 0.12
      else:
         argsn['bottom'] = 0.08
   if (argsn['right'] is None):
      if (ncol == 1):
         argsn['right'] = 0.85
      else:   
         argsn['right'] = 0.9
   if (argsn['left'] is None):
      if (ncol == 1):
         argsn['left'] = 0.12
      else:
         argsn['left'] = 0.08

   return argsn


def parse_plot_params(infname, ip, args):
   """
   Parse file in addition to command line
   """

   # create empty dictioary
   argsn = {}

   if (infname is None):
      input = ""
   else:
      # Read the input file into an array of strings -- one line for each element
      infile = open(infname,'r')
      input = infile.readlines()
      infile.close()

   # Input file must have syntax "parameter = value" with value being convertable
   # to appropriate type.
   
   argsn['fname'] = par_get_string('fname','panel{:d}'.format(ip),args['fname'],input)
   argsn['basename'] = par_get_string('basename','panel{:d}'.format(ip),args['basename'],input)
   argsn['outnum'] = par_get_int('outnum','panel{:d}'.format(ip),args['outnum'],input)
   argsn['var'] = par_get_string('var','panel{:d}'.format(ip),args['var'],input)
   argsn['label'] = par_get_string('label','panel{:d}'.format(ip),args['label'],input)
   argsn['cgs'] = par_get_boolean('cgs','panel{:d}'.format(ip),args['cgs'],input)
   argsn['type'] = par_get_string('type','panel{:d}'.format(ip),args['type'],input)

   # Parameters specific to pcolormesh plots
   argsn['colormap'] = par_get_string('colormap','panel{:d}'.format(ip),args['colormap'],input)
   argsn['ijplot'] = par_get_boolean('ijplot','panel{:d}'.format(ip),args['ijplot'],input)
   argsn['nocbar'] = par_get_boolean('nocbar','panel{:d}'.format(ip),args['nocbar'],input)
   argsn['logc'] = par_get_boolean('logc','panel{:d}'.format(ip),args['logc'],input)
   argsn['vlim'] = par_get_string('vlim','panel{:d}'.format(ip),args['vlim'],input)
   argsn['logr'] = par_get_string('logr','panel{:d}'.format(ip),args['logr'],input)

   # parameters for standard or both types of plots
   argsn['axis'] = par_get_int('axis','panel{:d}'.format(ip),args['axis'],input)
   argsn['i1'] = par_get_string('i1','panel{:d}'.format(ip),args['i1'],input)
   argsn['i2'] = par_get_string('i2','panel{:d}'.format(ip),args['i2'],input)
   argsn['xlim'] = par_get_string('xlim','panel{:d}'.format(ip),args['xlim'],input)
   argsn['ylim'] = par_get_string('ylim','panel{:d}'.format(ip),args['ylim'],input)
   argsn['logy'] = par_get_string('logy','panel{:d}'.format(ip),args['logy'],input)
   argsn['logx'] = par_get_string('logx','panel{:d}'.format(ip),args['logx'],input)

   return argsn

def limit_handler(limit):
   """
   Set axes limits
   """

   if (limit is None):
      return None
   elif (len(limit) > 1):
      # loop over all imu in the array
      slist = limit.strip(('[]')).split(",")
      lim = [float(i) for i in slist]

   if (len(lim) != 2):
      print("limit not formated correctly, returning None.")
      return None
   else:
      return lim

def index_handler(ind):
   """
   Parse index to deterimine which index to plot
   """

   if ind == None:
      return [0]
   elif ind == 'sum':
      return [ind]
   elif (len(ind) > 1):
      # loop over all indexes in the array
      slist = ind.strip(('[]')).split(",")
      ilist = [int(i) for i in slist]
   else:
      ilist = [ind]
   return ilist

def get_label(var,cgs,coord):
   """
   Determine label based on var
   """

   if (var == 'rho'):
      if (cgs):
         label = r'$\rho \; \rm{(g \; cm^{-3})}$'
      else:
         label = r'$\rho$'
   elif (var == 'tgas'): #MC variable
      if (cgs):
         label = r'$T \; \rm{(K)}$'
      else:
         label = r'$T$'
   elif (var == 'press'):
      if (cgs):
         label = r'$P \; \rm{(dyne)}$'
      else:
         label = r'$P$'
   elif (var == 'vel1'):
      if (coord == 'cartesian'):
         sub = r'_x'
      else:
         sub = r'_r'       
      if (cgs):
         label = r'$v'+sub+r' \; \rm{(cm \; s^{-1})}$'
      else:
         label = r'$v'+sub+r'$'
   elif (var == 'vel2'):
      if (coord == 'cartesian'):
         sub = r'_y'
      else:
         sub = r'_{\theta}'       
      if (cgs):
         label = r'$v'+sub+r' \; \rm{(cm \; s^{-1})}$'
      else:
         label = r'$v'+sub+r'$'
   elif (var == 'vel3'):
      if (coord == 'cartesian'):
         sub = r'_z'
      else:
         sub = r'_{\phi}'       
      if (cgs):
         label = r'$v'+sub+r' \; \rm{(cm \; s^{-1})}$'
      else:
         label = r'$v'+sub+r'$'
   else:
      label = None
   return label

def read_input_data(filename,basename=None,outnum=None,var=None,suffix='athdf',func=None):
   """
   Open the file to be read in and return dictionary
   """

   # set quantities
   if (var is None):
      if (func is None):
         print("Either var or func must be set. Exiting.")
         sys.exit(0)
      else:
         print("Fucntions not yet supported. Exiting.")
         sys.exit(0) 
   else:
      quantities = [var]

   # set input file
   if (filename is None):
      if (basename is None):
         print("Either filename or basename must be set. Exiting.")
         sys.exit(0)
      else:
         if (outnum is None):
            print("Basename set but outnum is not set. Exiting.")
            sys.exit(0)
         else:
            filename = filename = basename+".{:05d}".format(outnum)+"."+suffix
   else:
      if (basename is None):
         print("Both Basename and filename set. Defaulting to filename: "+filename)

   # read data into dictionary
   data = athena_read.athdf(filename,quantities=quantities)

   return data

def pcolormesh_one_panel(data,fig,gs,xaxis,yaxis,i3=0,**kwargs):
   """
   Plot one panel of figure using pcolormesh
   """
   
   #create new axis objects
   ax = fig.add_subplot(gs[0,0])
   axc = fig.add_subplot(gs[0,1]) #for colorbar

   # pop parameters
   #fname = kwargs.pop('fname')
   #basename = kwargs.pop('basename')
   #outnum = kwargs.pop('outnum')
   var = kwargs.pop('var')
   cgs = kwargs.pop('cgs')
   ijplot = kwargs.pop('ijplot')

   # Extract basic coordinate information
   coord = data['Coordinates'].decode('ascii', 'replace')
   if ((coord == 'spherical_polar') and (not ijplot)):
      r_f = data['x1f']
      if (kwargs['logr']):
         r_f = np.log(r_f)
      theta_f = data['x2f']
      
      # Make 2D mesh for plotting
      r_g, theta_g = np.meshgrid(r_f,theta_f)
      x_g = r_g * np.sin(theta_g)
      y_g = r_g * np.cos(theta_g)

      if (cgs):
         xaxis_label = r'$x \; \rm{(cm)}}$'
         yaxis_label = r'$y \; \rm{(cm)}}$'
      else:
         xaxis_label = r'$x$'
         yaxis_label = r'$y$'

   elif ((coord == 'spherical_polar') and ijplot):
      ic = np.arange(len(data['x1f']))
      jc = np.arange(len(data['x2f']))
      x_g, y_g = np.meshgrid(ic,jc)
      xaxis_label = r'$i$'
      yaxis_label = r'$j$'

   elif (coord == 'cartesian'):
      x = data['x1f']
      y = data['x2f']
      x_g, y_g = np.meshgrid(x,y)

      if (cgs):
         xaxis_label = r'$x \; \rm{(cm)}}$'
         yaxis_label = r'$y \; \rm{(cm)}}$'
      else:
         xaxis_label = r'$x$'
         yaxis_label = r'$y$'

   cmap = plt.get_cmap(kwargs['colormap'])

   # set vals assuming 2D for now
   vals = data[var][i3,:,:]
   if ijplot:
      vals = np.flip(vals,axis=0)
      
   if (kwargs['logc']):
      cnorm = colors.LogNorm()
   else:
      cnorm = colors.Normalize()
   
   vlim = limit_handler(kwargs['vlim'])
   if vlim is not None:
      vmin = vlim[0]
      vmax = vlim[1]
   else:
      vmin = None
      vmax = None

   pcm = ax.pcolormesh(x_g, y_g, vals, vmin=vmin, vmax=vmax, cmap=cmap, norm=cnorm)
   
   ax.set_xlabel(xaxis_label)
   ax.set_ylabel(yaxis_label)
   #ax.tick_params(labelsize=textsize)

   if (not xaxis):
      ax.set_xticklabels([])
   if (not yaxis):
      ax.set_yticklabels([])

   if (not kwargs['nocbar']):
      clabel = kwargs.pop('label')
      if (clabel is None):
         clabel = get_label(var,cgs,coord)
      cbar = Colorbar(ax = axc, mappable = pcm, orientation = 'vertical', ticklocation = 'right')
      cbar.set_label(clabel) 
      #cbar.ax.tick_params(labelsize=textsize) 

def plot_one_panel(data,fig,gs,xaxis,yaxis,**kwargs):
   """
   Plot one panel of standard plot
   """

   # create new axis object
   ax = fig.add_subplot(gs)

   # pop parameters
   #fname = kwargs.pop('fname')
   #basename = kwargs.pop('basename')
   #outnum = kwargs.pop('outnum')
   var = kwargs.pop('var')
   axis = kwargs.pop('axis')
   cgs = kwargs.pop('cgs')

   # Extract basic coordinate information
   coord = data['Coordinates'].decode('ascii', 'replace')

   # Reduce data to 1D curve and plot
   if (axis == 0):
      i2 = kwargs.pop('i2')
      ilist = index_handler(i2)
      #y = data[var][0,i2,:]
      x = data['x1v']
      for i in ilist:
         if (i == 'sum'):
            y = np.sum(data[var][0,:,:],axis=1)
         else:
            y = data[var][0,i,:]
         ax.plot(x,y)
      if (coord == 'cartesian'):
         xaxis_label = r'$x$'
      else:
         xaxis_label = r'$r$'
   elif (axis == 1):
      i1 = kwargs.pop('i1')
      ilist = index_handler(i1)
      #y = data[var][0,:,i1]
      x = data['x2v']
      for i in ilist:
         if (i == 'sum'):
            y = np.sum(data[var][0,:,:],axis=0)
         else:
            y = data[var][0,:,i]
         ax.plot(x,y)
      if (coord == 'cartesian'):
         xaxis_label = r'$y$'
      else:
         xaxis_label = r'$\theta$'

   # set axis scaling
   if (kwargs['logy']):
      ax.set_yscale('log')
   if (kwargs['logx']):
      ax.set_xscale('log')

   ylim = limit_handler(kwargs['ylim'])
   if ylim is not None:
      ax.set_ylim(ylim)
   xlim = limit_handler(kwargs['xlim'])
   if xlim is not None:
      ax.set_xlim(xlim)

   if (not xaxis):
      ax.set_xticklabels([])
   if (not yaxis):
      ax.set_yticklabels([])

   yaxis_label = kwargs.pop('label')
   if (yaxis_label is None):
      yaxis_label = get_label(var,cgs,coord)
   
   ax.plot(x,y)

   ax.set_xlabel(xaxis_label)
   ax.set_ylabel(yaxis_label) 
   #ax.tick_params(labelsize=textsize)

# Main function
def main(**kwargs):

   # set the figure geometry, overriding command line/defaults with content of pfile
   pfile = kwargs.pop('pfile')
   geo = set_plot_geometry(pfile,kwargs)
   npanels = geo['npanels']
   nrow = geo['nrow']
   ncol = geo['ncol']
   
   fig = plt.figure(figsize=(geo['figwidth'],geo['figheight']))
   rcParams.update({'font.size': geo['textsize']})

   width_ratios = [1] * ncol
   height_ratios = [1] * nrow

   grid = gridspec.GridSpec(nrow, ncol, width_ratios = width_ratios, height_ratios = height_ratios,
                            hspace = geo['hspace'], wspace = geo['wspace'], top = geo['top'], 
                            bottom = geo['bottom'], left = geo['left'], right = geo['right'])

   # loop over panels make nested gridspecs
   for i in range(npanels):
      j = i // ncol
      k = i % ncol

      # update arguments for this specific panel if requested in pfile
      args = parse_plot_params(pfile,i,kwargs)
      
      if (geo['axes_all']):
         xaxis = True
         yaxis = True
      else:
         if (j == nrow-1):
            xaxis = True
         else:
            xaxis = False
         if (k == 0):
            yaxis = True
         else:
            yaxis = False

      # read data and create arrays for pcolor mesh
      fname = args.pop('fname')
      basename = args.pop('basename')
      outnum = args.pop('outnum')
      data = read_input_data(fname,basename=basename,outnum=outnum,var=args['var'])

      # plot each panel
      if (args['type'] == 'pcolormesh'):
         # each panel is two subpanels, with one being the color bar
         gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec = grid[j,k], height_ratios = [1.],
                                               width_ratios = [1.,0.04] ,wspace = 0.05, hspace = 0.1)
         pcolormesh_one_panel(data,fig,gs,xaxis,yaxis,**args)
      elif (args['type'] == 'plot'):
         plot_one_panel(data,fig,grid[i],xaxis,yaxis,**args)

   plt.savefig("out.png")

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pfile',
                        default = None,
                        help='input file for additional parsing')
    parser.add_argument('-npanels',
                        default = 1,
                        type=int,
                        help='number of panels set via input file')
    parser.add_argument('-ncol',
                        type=int,
                        help='number of columns of panels')
    parser.add_argument('-axes_all',
                        action='store_true',
                        help='plot all panels with axes')
    parser.add_argument('-textsize',
                        default = None,
                        type=int,
                        help='textsize -- reset for many panel plots')
    parser.add_argument('-figheight',
                        default = None,
                        type=float,
                        help='height of figure in inches')
    parser.add_argument('-figwidth',
                        default = None,
                        type=float,
                        help='width of figure in inches')
    parser.add_argument('-hspace',
                        type=float,
                        default = None,
                        help='vertical whitespace bewtween panels')
    parser.add_argument('-wspace',
                        type=float,
                        default = None,
                        help='horizontal whitespace bewtween panel')
    parser.add_argument('-top',
                        type=float,
                        default = None,
                        help='maximum vertical extent')
    parser.add_argument('-bottom',
                        type=float,
                        default = None,
                        help='minimum vertical extent')
    parser.add_argument('-right',
                        type=float,
                        default = None,
                        help='maximum horizontal extent')
    parser.add_argument('-left',
                        type=float,
                        default = None,
                        help='minimum horizontal extent')
    # the following parameters may be set seperately for each panel
    parser.add_argument('-fname',
                        default=None,
                        help='input filename')
    parser.add_argument('-basename',
                        default=None,
                        help='input file basename')
    parser.add_argument('-outnum',
                        default = None,
                        type=int,
                        help='output number')
    parser.add_argument('-var',
                        default = None,
                        help='variable to be plot')
    parser.add_argument('-type',
                        default = 'pcolormesh',
                        help='type of plot')
    parser.add_argument('-label',
                        default = None,
                        help='label of plot variable')
    parser.add_argument('-cgs',
                        action='store_true',
                        help='set labels to cgs')
    parser.add_argument('-logc',
                        action='store_true',
                        help='flag indicating data should be colormapped logarithmically')
    parser.add_argument('-ijplot',
                        action='store_true',
                        help='treat r, theta as grids in i,j')
    parser.add_argument('-nocbar',
                        action='store_true',
                        help='plot without color bars')
    parser.add_argument('-c',
                        '--colormap',
                        default=None,
                        help=('name of Matplotlib colormap to use instead of default'))
    parser.add_argument('-axis',
                        default = 0,
                        type=int,
                        help='axis to use for dependent variable')
    parser.add_argument('-i1',
                        default = None,
                        help='slice # for x1 coordinate')
    parser.add_argument('-i2',
                        default = None,
                        help='slice # for x2 coordinate')
    parser.add_argument('-logr',
                        action='store_true',
                        help='use log scale for r axis in pcolormesh plots')
    parser.add_argument('-logx',
                        action='store_true',
                        help='use log scale for x axis')
    parser.add_argument('-logy',
                        action='store_true',
                        help='use log scale for y axis')
    parser.add_argument('-xlim',
                        default = None,
                        help='limits for x axis')
    parser.add_argument('-ylim',
                        default = None,
                        help='limits for y axis')
    parser.add_argument('-vlim',
                        default = None,
                        help='limits for color bar axis')
    args = parser.parse_args()
    main(**vars(args))
