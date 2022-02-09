#-------------------------------------------------------------------------------
# Module   : utils
# Purpose  : Miscellaneous low-level routines
# Author   : P. Rebischung
# Created  : 22-Aug-2014
#
# Changes  :
#
# Routines : - leap_sec      : Get number of leap seconds at requested dates
#            - temp_file     : Generate a random (UUID4) file name
#            - ps2png        : Convert PS file to PNG file
#            - sed_keywords  : Substitute keywords by their values in a string
#            - dict2rec      : Convert dictionary into record object
#            - rec2dict      : Convert record object into dictionary
#            - read_yaml     : Read YAML configuration file
#            - write_yaml    : Write YAML configuration file
#            - parallel_sh   : Execute multiple commands in parallel
#            - print_message : Print a message both to screen and to log file
#-------------------------------------------------------------------------------

# LIBRARIES
#-------------------------------------------------------------------------------
import os
import re
import uuid
import yaml
import numpy

from .constants import *
from .date import date


#-------------------------------------------------------------------------------
# Routine : leapsec
# Purpose : Get number of leap seconds at requested dates
# Author  : P. Rebischung
# Created : 16-Dec-2011
#
# Changes :
#
# Input   : - d         : Requested dates (MJD)
#           - timescale : 'GPS' or 'UTC' depending on time scale of d. Default is 'UTC'.
# Output  : - l         : GPS-UTC leap seconds at requested dates
#-------------------------------------------------------------------------------
def leapsec(d, timescale='UTC'):

  # Initialization
  l = numpy.zeros(len(d))

  # 1st case : MJDs given in GPS time
  if (timescale == 'GPS'):
    for i in range(len(d)):
      ind = numpy.nonzero(mjd_leap <= d[i])[0][-1]
      l[i] = gps_utc[ind]

  # 2nd case : MJDs given in UTC
  elif (timescale == 'UTC'):
    for i in range(len(d)):
      ind = numpy.nonzero(mjd_leap - gps_utc / 86400 <= d[i])[0][-1]
      l[i] = gps_utc[ind]

  return l


#-------------------------------------------------------------------------------
# Routine : temp_file
# Purpose : Generate a random (UUID4) file name
# Author  : P. Rebischung
# Created : 27-Aug-2014
#
# Changes :
#
# Input   : 
# Output  : 
#-------------------------------------------------------------------------------
def temp_file():
  
  return str(uuid.uuid4())


#-------------------------------------------------------------------------------
# Routine : ps2png
# Purpose : Convert PS file to PNG file
# Author  : P. Rebischung
# Created : 28-Feb-2012
#
# Changes :
#
# Input   : - ps     : PS file
#           - png    : PNG file
#           - rotate : Rotation angle. Default is 0.
#           - margin : Margin. Default is 20.
# Output  : 
#-------------------------------------------------------------------------------
def ps2png(ps, png, rotate=0, margin=20):

  # Temporary file
  tmp = temp_file() + '.png'

  # Let's go!
  os.system('gs -dQUIET -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r250 -dGraphicsAlphaBits=4 -sOutputFile={0} {1}'.format(tmp, ps))
  os.system('convert {0} -rotate {1} -quality 100 -trim -mattecolor white -frame {2}x{2} {3}'.format(tmp, rotate, margin, png))
  os.system('rm {0}'.format(tmp))

  
#-------------------------------------------------------------------------------
# Routine : sed_keywords
# Purpose : Substitute keywords by their values in a string
# Author  : P. Rebischung
# Created : 22-Aug-2014
#
# Changes :
#
# Input   : - s : String
#           - t : date object
# Output  : - s : Updated string
#-------------------------------------------------------------------------------
def sed_keywords(s, t):
  
  # Date elements
  s = re.subn('\$yyyy', t.yyyy, s)[0]
  s = re.subn('\$doy' , t.doy, s)[0]
  s = re.subn('\$yy', t.yy, s)[0]
  s = re.subn('\$mm', t.mm, s)[0]
  s = re.subn('\$dd', t.dd, s)[0]
  s = re.subn('\$hour', t.hour, s)[0]
  s = re.subn('\$min' , t.min, s)[0]
  s = re.subn('\$sec' , t.sec, s)[0]
  s = re.subn('\$week', t.week, s)[0]
  s = re.subn('\$dow' , t.dow, s)[0]
  s = re.subn('\$wk', t.wk, s)[0]
  
  # Operating system
  s = re.subn('\$os', os.uname()[0], s)[0]
  
  return s


#-------------------------------------------------------------------------------
# Routine : dict2rec
# Purpose : Convert dictionary into record object
# Author  : P. Rebischung
# Created : 03-Mar-2015
#
# Changes : 
#
# Input   : - y   : Dictionary
#           - sed : If True, make use of sed_keywords for each key. Default is False.
#           - t   : date object required if sed=True. Default is None.
# Output  : - r   : record object
#-------------------------------------------------------------------------------
def dict2rec(y, sed=False, t=None):
  
  # Initialization
  r = record()
  
  # Loop over dictionary keys
  for key in y:
    
    # If y[key] is a dictionary,
    if (type(y[key]).__name__ == 'dict'):
      setattr(r, key, dict2rec(y[key], sed, t))
      
    # If y[key] is an empty list,
    elif ((type(y[key]).__name__ == 'list') and (len(y[key]) == 0)):
      setattr(r, key, y[key])
      
    # If y[key] is a list of dictionaries,
    elif ((type(y[key]).__name__ == 'list') and (type(y[key][0]).__name__ == 'dict')):
      setattr(r, key, [dict2rec(d, sed, t) for d in y[key]])
      
    # If y[key] is a list of strings and sed is needed,
    elif ((type(y[key]).__name__ == 'list') and (type(y[key][0]).__name__ == 'str') and (sed)):
      setattr(r, key, [sed_keywords(s, t) for s in y[key]])
      
    # If y[key] is a string and sed is needed,
    elif ((type(y[key]).__name__ == 'str') and (sed)):
      setattr(r, key, sed_keywords(y[key], t))
    
    # Other cases
    else:
      setattr(r, key, y[key])
      
  return r


#-------------------------------------------------------------------------------
# Routine : rec2dict
# Purpose : Convert record object into dictionary
# Author  : P. Rebischung
# Created : 03-Mar-2015
#
# Changes : 
#
# Input   : - r : record object
# Output  : - d : dictionary
#-------------------------------------------------------------------------------
def rec2dict(r):
  
  # Initialization
  d = vars(r)
  
  # Loop over dictionary keys
  for key in d:
    
    # If d[key] is a record,
    if (type(d[key]).__name__ == 'instance'):
      d[key] = rec2dict(d[key])
      
    # If d[key] is a non-empty list of records,
    elif ((type(d[key]).__name__ == 'list') and (len(d[key]) > 0)):
      if (type(d[key][0]).__name__ == 'instance'):
        d[key] = [rec2dict(r) for r in d[key]]
        
  return d

  
#-------------------------------------------------------------------------------
# Routine : read_yaml
# Purpose : Read YAML configuration file
# Author  : P. Rebischung
# Created : 10-Oct-2013
#
# Changes : PR, 03-Mar-2015 : Make recursive use of subroutine dict2rec
#
# Input   : - inp   : YAML configuration file
#           - sed   : If True, make use of sed_keywords for each key. Default is False.
#           - t     : date object required if sed=True. Default is None.
# Output  : Record or list
#-------------------------------------------------------------------------------
def read_yaml(inp, sed=False, t=None):
  
    y = yaml.load(open(inp, 'r'))
    
    # If y is a dictionary,
    if (type(y).__name__ == 'dict'):
      return dict2rec(y, sed, t)
    
    # If y is an empty list,
    elif ((type(y).__name__ == 'list') and (len(y) == 0)):
      return y
    
    # If y is a list of dictionaries,
    elif ((type(y).__name__ == 'list') and (type(y[0]).__name__ == 'dict')):
      return [dict2rec(d, sed, t) for d in y]
        
    # If y is a list of strings and sed is needed,
    elif ((type(y).__name__ == 'list') and (type(y[0]).__name__ == 'str') and (sed)):
      return [sed_keywords(s, t) for s in y]
    
    # Other cases
    else:
      return y


#-------------------------------------------------------------------------------
# Routine : write_yaml
# Purpose : Write YAML configuration file
# Author  : P. Rebischung
# Created : 30-Sep-2014
#
# Changes : PR, 03-Mar-2015 : Make recursive use of subroutine rec2dict
#
# Input   : - r   : Record or list
#           - out : File to write
# Output  : 
#-------------------------------------------------------------------------------
def write_yaml(r, out):
  
  # If r is a record,
  if (type(r).__name__ == 'instance'):
    l = rec2dict(r)
    
  # If r is an empty list,
  elif ((type(r).__name__ == 'list') and (len(r) == 0)):
    l = r
  
  # If r is a list of records,
  elif ((type(r).__name__ == 'list') and (type(r[0]).__name__ == 'instance')):
    l = [rec2dict(rec) for rec in r]
      
  # Other cases
  else:
    l = r
      
  yaml.dump(l, (out, 'w'))


#-------------------------------------------------------------------------------
# Routine : parallel_sh
# Purpose : Execute multiple commands in parallel
# Author  : P. Rebischung
# Created : 25-Aug-2014
#
# Changes :
#
# Input   : - file  : File containing the list of commands to execute
#           - nproc : Number of processors to use
# Output  : 
#-------------------------------------------------------------------------------
def parallel_sh(file, nproc):
  
  # Read list of commands
  commands = open(file).readlines()
  
  # Write make file
  makefile = temp_file() + '.make'
  f = open(makefile, 'w')
  f.write('all :')
  for i in range(len(commands)):
    f.write(' job{0}'.format(i))
  f.write('\n')
  for i in range(len(commands)):
    f.write('job{0} :\n\t{1}'.format(i, commands[i]))
  f.close()
  
  # Execute make file
  os.system('make -j{0} -f {1}'.format(nproc, makefile))
  
  # Remove make file
  os.system('rm {0}'.format(makefile))


#-------------------------------------------------------------------------------
# Routine : print_message
# Purpose : Print a message both to screen and to log file
# Author  : P. Rebischung
# Created : 12-May-2012
#
# Changes :
#
# Input   : - message   : Message
#           - mainlog   : Log file handle
#           - timestamp : True if a time stamp should be added at the beginning of message
#-------------------------------------------------------------------------------
def print_message(message, mainlog, timestamp=False):

  # Add time stamp if requested
  if (timestamp):
    message = str(date()) + ' : ' + message

  # Print message to screen
  print(message)

  # Print message to main log file
  mainlog.write(message + '\n')

  
#-------------------------------------------------------------------------------
# Routine : isfloat
# Purpose : Test if a string can be converted into float
# Author  : P. Rebischung
# Created : 02-Feb-2017
#
# Changes :
#
# Input   : s : String
# Output  : True or False
#-------------------------------------------------------------------------------
def isfloat(s):
  
  try:
    float(s)
    return True
  
  except ValueError:
    return False
