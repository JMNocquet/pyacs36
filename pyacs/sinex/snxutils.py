#-------------------------------------------------------------------------------
# Module   : snxutils
# Purpose  : Low-level routines for SINEX manipulation
# Author   : P. Rebischung
# Created  : 22-Aug-2014
#
# Changes  :
#
# Routines : - earlier            : Compare two dates in SINEX format
#            - read_domes         : Read DOMES number catalogue (codomes.snx)
#            - read_solns         : Read discontinuity list (soln.snx)
#            - extract_solns      : Extract a set of stations from a discontinuity list
#            - helmert_outliers   : Find outliers in the residuals of a 7-parameter Helmert comparison
#            - helmert14_outliers : Find outliers in the residuals of a 14-parameter Helmert comparison
#            - helmert_hist       : Draw histograms of Helmert residuals
#            - helmert_map        : Draw map of Helmert residuals
#            - station_map        : Draw map of stations
#            - meanpole           : Compute long term mean pole
#            - read_opoleload     : Read grid of ocean pole tide loading coefficients
#            - compute_opoleload  : Read grid of ocean pole tide loading coefficients
#            - compute_psd        : Compute post-seismic deformations and associated uncertainties
#-------------------------------------------------------------------------------



# LIBRARIES
#-------------------------------------------------------------------------------
import os
import matplotlib.pyplot as pp
import matplotlib.ticker as ticker
#from mpl_toolkits.basemap import Basemap

from .constants import *
from .mathutils import *
from .date import date



#-------------------------------------------------------------------------------
# Routine : earlier
# Purpose : Compare two dates in SINEX format
# Author  : P. Rebischung
# Created : 12-Aug-2011
#
# Changes :
#
# Input   : - t1 : Date in sinex format
#           - t2 : Date in sinex format
# Output  : - b  : True if t1 <= t2. False otherwise.
#-------------------------------------------------------------------------------
def earlier(t1, t2):

  if (int(t1[0:2]) >= 50):
    t1 = '19' + t1
  else:
    t1 = '20' + t1

  if (int(t2[0:2]) >= 50):
    t2 = '19' + t2
  else:
    t2 = '20' + t2

  return (t1 <= t2)

#-------------------------------------------------------------------------------
# Routine : read_domes
# Purpose : Read DOMES number catalogue (codomes.snx)
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes : PR, 27-Jan-2017 : Add possibility to read station coordinates
#
# Input   : - file  : DOMES number catalogue
#           - coord : True if station coordinates should be read. Default is True.
# Output  : - domes : List of records
#-------------------------------------------------------------------------------
def read_domes(file, coord=True):

  # Initialization
  domes = []

  # Open input file
  f = open(file, 'rb')

  # Read file
  line = f.readline()
  while (line):
    r = record()
    r.code = line[0:4]
    r.pt = line[5:7]
    r.domes = line[8:17]
    r.description = line[19:41]
    
    if (coord):
      tab = line.strip().split()
      r.lon = float(tab[-2])
      r.lat = float(tab[-1])
    
    domes.append(r)
    line = f.readline()

  # Some manual changes to handle duplicate DOMES numbers
  for i in range(len(domes)):
    if (domes[i].code == 'GOLD'):
      domes[i].domes = '40405S031'
    elif (domes[i].code == 'IISC'):
      domes[i].domes = '22306M002'
    elif (domes[i].code == 'KELY'):
      domes[i].domes = '43005M002'
    elif (domes[i].code == 'MDVO'):
      domes[i].domes = '12309M002'
    elif (domes[i].code == 'MTKA'):
      domes[i].domes = '21741S002'
    elif (domes[i].code == 'UPAD'):
      domes[i].domes = '12750M002'
    elif (domes[i].code == 'WEL2'):
      domes[i].domes = '50208S003'

  # Close file
  f.close()

  return domes

#-------------------------------------------------------------------------------
# Routine : read_solns
# Purpose : Read discontinuity list (soln.snx)
# Author  : P. Rebischung
# Created : 11-Aug-2011
#
# Changes :
#
# Input   : - file  : Discontinuity list
# Output  : - solns : List of records
#-------------------------------------------------------------------------------
def read_solns(file):

  # Initialization
  solns = []

  code  = ''
  pt    = ''
  ista  = -1    

  # Open input file
  f = open(file, 'rb')

  # Skip header
  line = f.readline()
  while (line[0:23] != '+SOLUTION/DISCONTINUITY'):
    line = f.readline()

  # Read SOLUTION/DISCONTINUITY block
  line = f.readline()
  while (line[0:1] != '-'):
    if (line[0:1] != '*'):

      # New station?
      if ((line[1:5] != code) or (line[6:8] != pt)):
        ista = ista+1
        code = line[1:5]
        pt   = line[6:8]

        r = record()
        r.code = code
        r.pt   = pt
        r.P    = []
        r.V    = []
        solns.append(r)

      # Position soln?
      if (line[42:43] == 'P'):
        r = record()
        r.soln   = line[9:13]
        r.start  = line[16:28]
        r.end    = line[29:41]
        r.reason = line[46:].strip()
        solns[ista].P.append(r)

      # Velocity soln?
      elif (line[42:43] == 'V'):
        r = record()
        r.soln   = line[9:13]
        r.start  = line[16:28]
        r.end    = line[29:41]
        r.reason = line[46:].strip()
        solns[ista].V.append(r)

    line = f.readline()

  # Close file
  f.close()

  return solns

#-------------------------------------------------------------------------------
# Routine : extract_solns
# Purpose : Extract a set of stations from a discontinuity list
# Author  : P. Rebischung
# Created : 25-Jun-2012
#
# Changes : 04-Oct-2012, PR : Small bugs corrected
#
# Input   : - solns   : Discontinuity list
#           - snx     : sinex object defining the list of stations to extract
#           - file    : Output discontinuity file
#           - comment : File containing FILE/COMMENT block for the output file. Default if None.
# Output  : 
#-------------------------------------------------------------------------------
def extract_solns(solns, snx, file, comment=None):

  # List of stations to be extracted
  codept = [s.code+s.pt for s in snx.sta]

  # New solns table
  solns2 = []
  for s in solns:
    if (s.code+s.pt in codept):
      solns2.append(s)

  # Write 1st line of output file
  f = open(file, 'w')
  f.write('%=SNX 2.01 {0} {1} {0} 00:000:00000 00:000:00000 P     0 1 X V\n'.format(agency, date().snxepoch()))
  f.close()

  # Cat FILE/COMMENT block if any
  if (comment):
    os.system('cat {0} >> {1}'.format(comment, file))

  # Write SOLUTION/DISCONTINUITY block
  f = open(file, 'a')
  f.write('+SOLUTION/DISCONTINUITY\n')
  
  for s in solns2:
    for p in s.P:
      f.write(' {0.code} {0.pt} {1.soln} P {1.start} {1.end} P - {1.reason}\n'.format(s, p))
    for v in s.V:
      f.write(' {0.code} {0.pt} {1.soln} P {1.start} {1.end} V - {1.reason}\n'.format(s, v))
      
  f.write('-SOLUTION/DISCONTINUITY\n')
  f.write('*----------------------------------------------------------------\n')
  f.write('%ENDSNX\n')
  f.close()

#-------------------------------------------------------------------------------
# Routine : helmert_outliers
# Purpose : Find outliers in the residuals of a 7-parameter Helmert comparison
# Author  : P. Rebischung
# Created : 07-Jul-2011
#
# Changes :
#
# Input   : - code      : 4-char codes
#           - pt        : PT codes
#           - soln      : Solns
#           - v         : Raw residuals (mm)
#           - vn        : Normalized residuals
#           - thr_raw   : Table of [E, N, H] thresholds in mm for raw residuals. Default is None.
#           - thr_norm  : Table of [E, N, H] thresholds for normalized residuals. Default is None.
#           - thr_auto  : If neither thr_raw nor thr_norm are specified, the "automatic" thresholds will be
#                         thr_auto*[RMS_E, RMS_N, RMS_H] for both raw and normalized residuals. Default is 4.
#           - log       : Log file. Default is None.
#           - append    : If true, text will be appended to log file. Default if False.
# Output  : - codedel   : 4-char codes of outliers
#           - ptdel     : PT codes of outliers
#           - solndel   : Solns of outliers
#-------------------------------------------------------------------------------
def helmert_outliers(code, pt, soln, v, vn, thr_raw=None, thr_norm=None, thr_auto=4., log=None, append=False):

  # Define thresholds in case of an automatic detection
  if ((thr_raw == None) and (thr_norm == None)):
    thr_raw  = [thr_auto*rms(v[:,0]),  thr_auto*rms(v[:,1]), thr_auto*rms(v[:,2])]
    thr_norm = [thr_auto*rms(vn[:,0]), thr_auto*rms(vn[:,1]), thr_auto*rms(vn[:,2])]
    
  # Get indices of outliers
  ind_raw  = numpy.array([], dtype='i8')
  ind_norm = numpy.array([], dtype='i8')
  if (thr_raw != None):
    ind_raw = numpy.nonzero(numpy.any(numpy.abs(v) > thr_raw, axis=1))[0]
  if (thr_norm != None):
    ind_norm = numpy.nonzero(numpy.any(numpy.abs(vn) > thr_norm, axis=1))[0]
  ind = numpy.union1d(ind_raw, ind_norm)

  # Get outlier identifiers
  codedel = [code[i] for i in ind]
  ptdel   = [pt[i]   for i in ind]
  solndel = [soln[i] for i in ind]

  # Write log file if necessary
  if (log):
    if (append):
      f = open(log, 'a')
    else:
      f = open(log, 'w')

    # Write header
    f.write('================================================================================\n')
    f.write('snxutils::helmert_outliers : Find outliers in 7-parameter Helmert residuals\n')
    f.write('================================================================================\n')
    f.write('\n')

    # Print used thresholds
    f.write('Used thresholds :\n')
    f.write('-----------------\n')
    f.write('\n')
    f.write('  Raw thresholds (mm) | Normalized threholds |\n')
    f.write('----------------------|----------------------|\n')
    f.write('    E      N      H   |    E      N      H   |\n')
    f.write('----------------------|----------------------|\n')
    if (thr_raw != None):
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |'.format(thr_raw))
    else:
      f.write('                      |')
    if (thr_norm != None):
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |\n'.format(thr_norm))
    else:
      f.write('                      |\n')
    f.write('----------------------|----------------------|\n')
    f.write('\n')

    # Print detected outliers
    f.write('Detected outliers :\n')
    f.write('-------------------\n')
    f.write('\n')
    f.write('             |  Raw residuals (mm)  | Normalized residuals |\n')
    f.write('-------------|----------------------|----------------------|\n')
    f.write('code pt soln |    E      N      H   |    E      N      H   |\n')
    f.write('-------------|----------------------|----------------------|\n')
    for i in range(len(ind)):
      f.write('{0} {1} {2} |'.format(codedel[i], ptdel[i], solndel[i]))
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |'.format(v[ind[i]]))
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} |\n'.format(vn[ind[i]]))
    f.write('-------------|----------------------|----------------------|\n')
    f.write('\n')

    f.close()

  return (codedel, ptdel, solndel)

#-------------------------------------------------------------------------------
# Routine : helmert14_outliers
# Purpose : Find outliers in the residuals of a 14-parameter Helmert comparison
# Author  : P. Rebischung
# Created : 07-Jul-2011
#
# Changes :
#
# Input   : - code      : 4-char codes
#           - pt        : PT codes
#           - soln      : Solns
#           - v         : Raw residuals (mm, mm/y)
#           - vn        : Normalized residuals
#           - thr_raw   : Table of [E, N, H, VE, VN, VH] thresholds in mm[/y] for raw residuals. Default is None.
#           - thr_norm  : Table of [E, N, H, VE, VN, VH] thresholds for normalized residuals. Default is None.
#           - thr_auto  : If neither thr_raw nor thr_norm are specified, the "automatic" thresholds will be
#                         thr_auto*[RMS_E, RMS_N, RMS_H, RMS_VE, RMS_VN, RMS_VH] for both raw and normalized
#                         residuals. Default is 5.
#           - log       : Log file. Default is None.
#           - append    : If true, text will be appended to log file. Default if False.
# Output  : - codedel   : 4-char codes of outliers
#           - ptdel     : PT codes of outliers
#           - solndel   : Solns of outliers
#-------------------------------------------------------------------------------
def helmert14_outliers(code, pt, soln, v, vn, thr_raw=None, thr_norm=None, thr_auto=5.0, log=None, append=False):

  # Define thresholds in case of an automatic detection
  if ((not(thr_raw)) and (not(thr_norm))):
    thr_raw  = [thr_auto*rms(v[:,0]),  thr_auto*rms(v[:,1]), thr_auto*rms(v[:,2]), thr_auto*rms(v[:,3]),  thr_auto*rms(v[:,4]), thr_auto*rms(v[:,5])]
    thr_norm = [thr_auto*rms(vn[:,0]), thr_auto*rms(vn[:,1]), thr_auto*rms(vn[:,2]), thr_auto*rms(vn[:,3]), thr_auto*rms(vn[:,4]), thr_auto*rms(vn[:,5])]
    
  # Get indices of outliers
  ind_raw  = numpy.array([], dtype='i8')
  ind_norm = numpy.array([], dtype='i8')
  if (thr_raw):
    ind_raw = numpy.nonzero(numpy.any(numpy.abs(v) > thr_raw, axis=1))[0]
  if (thr_norm):
    ind_norm = numpy.nonzero(numpy.any(numpy.abs(vn) > thr_norm, axis=1))[0]
  ind = numpy.union1d(ind_raw, ind_norm)

  # Get outlier identifiers
  codedel = [code[i] for i in ind]
  ptdel   = [pt[i]   for i in ind]
  solndel = [soln[i] for i in ind]

  # Write log file if necessary
  if (log):
    if (append):
      f = open(log, 'a')
    else:
      f = open(log, 'w')

    # Write header
    f.write('================================================================================\n')
    f.write('snxutils::helmert14_outliers : Find outliers in 14-parameter Helmert residuals\n')
    f.write('================================================================================\n')
    f.write('\n')

    # Print used thresholds
    f.write('Used thresholds :\n')
    f.write('-----------------\n')
    f.write('\n')              
    f.write('         Raw thresholds (mm, mm/y)         |            Normalized threholds           |\n')
    f.write('-------------------------------------------|-------------------------------------------|\n')
    f.write('    E      N      H     VE     VN     VH   |    E      N      H     VE     VN     VH   |\n')
    f.write('-------------------------------------------|-------------------------------------------|\n')
    if (thr_raw):
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} {0[3]:6.1f} {0[4]:6.1f} {0[5]:6.1f} |'.format(thr_raw))
    else:
      f.write('                                           |')
    if (thr_norm):
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} {0[3]:6.1f} {0[4]:6.1f} {0[5]:6.1f} |\n'.format(thr_norm))
    else:
      f.write('                                           |\n')
    f.write('-------------------------------------------|-------------------------------------------|\n')
    f.write('\n')

    # Print detected outliers
    f.write('Detected outliers :\n')
    f.write('-------------------\n')
    f.write('\n')
    f.write('             |         Raw residuals (mm, mm/y)          |           Normalized residuals            |\n')
    f.write('-------------|-------------------------------------------|-------------------------------------------|\n')
    f.write('code pt soln |    E      N      H     VE     VN     VH   |    E      N      H     VE     VN     VH   |\n')
    f.write('-------------|-------------------------------------------|-------------------------------------------|\n')
    for i in range(len(ind)):
      f.write('{0} {1} {2} |'.format(codedel[i], ptdel[i], solndel[i]))
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} {0[3]:6.1f} {0[4]:6.1f} {0[5]:6.1f} |'.format(v[ind[i]]))
      f.write(' {0[0]:6.1f} {0[1]:6.1f} {0[2]:6.1f} {0[3]:6.1f} {0[4]:6.1f} {0[5]:6.1f} |\n'.format(vn[ind[i]]))
    f.write('-------------|-------------------------------------------|-------------------------------------------|\n')
    f.write('\n')

    f.close()

  return(codedel, ptdel, solndel)

#-------------------------------------------------------------------------------
# Routine : helmert_hist
# Purpose : Draw histograms of Helmert residuals
# Author  : P. Rebischung
# Created : 08-Aug-2011
#
# Changes : PR, 13-Jun-2016: Add option 'no_up'
#
# Input   : - v      : Residuals
#           - unit   : Residuals unit. Default is 'mm'.
#           - nbins  : Number of bins for each histogram. Default is 100.
#           - output : Output file. Default is None.
#           - no_up  : True if only horizontal residuals. Default is False.
# Output  : 
#-------------------------------------------------------------------------------
def helmert_hist(v, unit='mm', nbins=100, output=None, no_up=False):

  # Set X-axis limit
  xmax = ceil(numpy.max(numpy.abs(v)))

  # Bins
  width = 2.*xmax/nbins
  bins = numpy.arange(-xmax, xmax+width, width)
  
  # Set Y-axis limit
  (he, bins) = numpy.histogram(v[:,0], bins)
  (hn, bins) = numpy.histogram(v[:,1], bins)
  if (not(no_up)):
    (hu, bins) = numpy.histogram(v[:,2], bins)
    ymax = 1.1 * numpy.max([he, hn, hu])
  else:
    ymax = 1.1 * numpy.max([he, hn])

  # Set X-axis annotation spacing
  if (xmax <= 2):
    ax = 0.5
  elif (xmax <= 5):
    ax = 1
  elif (xmax <= 10):
    ax = 2.5
  elif (xmax <= 20):
    ax = 5
  elif (xmax <= 50):
    ax = 10
  elif (xmax <= 100):
    ax = 25
  elif (xmax <= 200):
    ax = 50
  elif (xmax <= 500):
    ax = 100
  elif (xmax <= 1000):
    ax = 250
  elif (xmax <= 2000):
    ax = 500
  else:
    ax = 1000

  # Set Y-axis annotation spacing
  if (ymax <= 5):
    ay = 1
  elif (ymax <= 10):
    ay = 2
  elif (ymax <= 25):
    ay = 5
  elif (ymax <= 50):
    ay = 10
  elif (ymax <= 100):
    ay = 20
  elif (ymax <= 250):
    ay = 50
  elif (ymax <= 500):
    ay = 100
  elif (ymax <= 1000):
    ay = 200
  elif (ymax <= 2500):
    ay = 500
  else:
    ay = 1000

  # Create new figure and define tick options
  fig = pp.figure(figsize=(6, 12))
  formatter = ticker.FormatStrFormatter('%.0f {0}'.format(unit))
  xlocator  = ticker.MultipleLocator(ax)
  ylocator  = ticker.MultipleLocator(ay)

  # Draw histogram of East residuals
  sub1 = fig.add_subplot(311)
  sub1.hist(v[:,0], bins)
  sub1.set_ylabel('# stations')
  sub1.axis([-xmax, xmax, 0, ymax])
  sub1.xaxis.set_major_locator(xlocator)
  sub1.yaxis.set_major_locator(ylocator)
  sub1.xaxis.set_major_formatter(formatter)
  sub1.set_title('East residuals', position=(0.5, 1.06))

  # Draw histogram of North residuals
  sub2 = fig.add_subplot(312)
  sub2.hist(v[:,1], bins)
  sub2.set_ylabel('# stations')
  sub2.axis([-xmax, xmax, 0, ymax])
  sub2.xaxis.set_major_locator(xlocator)
  sub2.yaxis.set_major_locator(ylocator)
  sub2.xaxis.set_major_formatter(formatter)
  sub2.set_title('North residuals', position=(0.5, 1.06))

  # Draw histogram of Up residuals
  if (not(no_up)):
    sub3 = fig.add_subplot(313)
    sub3.hist(v[:,2], bins)
    sub3.set_ylabel('# stations')
    sub3.axis([-xmax, xmax, 0, ymax])
    sub3.xaxis.set_major_locator(xlocator)
    sub3.yaxis.set_major_locator(ylocator)
    sub3.xaxis.set_major_formatter(formatter)
    sub3.set_title('Up residuals', position=(0.5, 1.06))
  
  # Compute ENH RMSs
  rmse = rms(v[:,0])
  rmsn = rms(v[:,1])
  if (not(no_up)):
    rmsu = rms(v[:,2])

  # Print RMSs
  sub1.text(0.02, 0.95, 'RMS = {0:.1f} {1}'.format(rmse, unit), transform=sub1.transAxes, va='top', ha='left')
  sub2.text(0.02, 0.95, 'RMS = {0:.1f} {1}'.format(rmsn, unit), transform=sub2.transAxes, va='top', ha='left')
  if (not(no_up)):
    sub3.text(0.02, 0.95, 'RMS = {0:.1f} {1}'.format(rmsu, unit), transform=sub3.transAxes, va='top', ha='left')

  # Fit and draw gaussians
  ae = len(v) * width / (sqrt(2*pi) * rmse)
  an = len(v) * width / (sqrt(2*pi) * rmsn)
  if (not(no_up)):
    au = len(v) * width / (sqrt(2*pi) * rmsu)

  x = numpy.arange(-xmax, 1.002*xmax, 0.002*xmax)
  ye = ae * numpy.exp(-x**2 / (2*rmse**2))
  yn = an * numpy.exp(-x**2 / (2*rmsn**2))
  if (not(no_up)):
    yu = au * numpy.exp(-x**2 / (2*rmsu**2))

  sub1.plot(x, ye, 'r', linewidth=2)
  sub2.plot(x, yn, 'r', linewidth=2)
  if (not(no_up)):
    sub3.plot(x, yu, 'r', linewidth=2)

  # Save figure into output file...
  if (output):
    pp.savefig(output, bbox_inches='tight')
    pp.close()

  # ...or show it.
  else:
    pp.show()

#-------------------------------------------------------------------------------
# Routine : helmert_map
# Purpose : Draw residuals of a 7-parameter Helmert comparison on a map
# Author  : P. Rebischung
# Created : 09-Aug-2011
#
# Changes :
#
# Input   : - lon     : Longitudes
#           - lat     : Latitudes
#           - v       : Residuals
#           - unit    : Residuals unit. Default is 'mm'.
#           - h_scale : Scale factor for horizontal residuals (wrt default scale).
#           - v_scale : Scale factor for vertical residuals (wrt default scale).
#           - legend  : Length of legend vectors in mm. Default is 5.
#           - title   : Map title. Default is None.
#           - output  : Output file. Default is None.
#           - no_up   : True if only horizontal residuals. Default is False.
# Output  : 
#-------------------------------------------------------------------------------
def helmert_map(lon, lat, v, unit='mm', h_scale=1, v_scale=1, legend=5, title=None,
                output=None, no_up=False):

  # Append legend points to lon, lat and v lists
  if (not(no_up)):
    lon = numpy.hstack((lon, [-140, -141, -141]))
    lat = numpy.hstack((lat, [-42, -52, -52]))
    v = numpy.vstack((v, [[legend, 0, 0], [0, 0, legend], [0, 0, -legend]]))
  else:
    lon = numpy.hstack((lon, -140))
    lat = numpy.hstack((lat, -42))
    v = numpy.vstack((v, [legend, 0]))

  # Draw basemap
  pp.figure(figsize=(16, 12))
  map = Basemap(projection='robin', lon_0=0, resolution='c', area_thresh=50000.)
  map.drawmapboundary(fill_color=(0.82, 0.86, 1, 1)) 
  map.fillcontinents(color='white', lake_color=(0.82, 0.86, 1, 1), zorder=1)
  map.drawcoastlines(color='blue', linewidth=0.3)

  # Add title if necessary
  if (title):
    pp.title(title, position=(0.5, 1.06))

  # Plot points
  (x, y) = list(map(lon, lat))
  map.plot(x, y, '.', color='black', markersize=6)

  # Rotate horizontal residuals into map projection
  lon2 = lon + 180*v[:,0] / (pi*ae * numpy.cos(pi/180*lat))
  lat2 = lat + 180*v[:,1] / (pi*ae)
  (x2, y2) = list(map(lon2, lat2))
  vx = x2 - x
  vy = y2 - y

  # Plot horizontal residuals
  map.quiver(x, y, vx, vy, zorder=2, units='dots', width=1.5, scale=0.1/h_scale, color='black')

  # Plot vertical residuals
  if (not(no_up)):
    for i in range(len(lon)):
      if (v[i,2] > 0):
        map.plot([x[i], x[i]], [y[i], y[i] + 1.5e5*v_scale*v[i,2]], color='red', linewidth=1.5)
      else:
        map.plot([x[i], x[i]], [y[i], y[i] + 1.5e5*v_scale*v[i,2]], color='green', linewidth=1.5)

  # Legend text for horizontal residuals
  (xt, yt) = list(map(-101, -42))
  pp.text(xt, yt, '{0} {1}'.format(legend, unit), ha='right', va='center')

  # Legend text for vertical residuals
  if (not(no_up)):
    (xt, yt) = list(map(-107, -52))
    pp.text(xt, yt, r'$\pm$ {0} {1}'.format(legend, unit), ha='right', va='center')

  # Legend box
  if (not(no_up)):
    (xb, yb) = list(map([-142, -95, -113, -169, -142], [-38, -38, -62, -62, -38]))
    map.plot(xb, yb, 'k')
  else:
    (xb, yb) = list(map([-142, -95, -99.2, -148.4, -142], [-38, -38, -46, -46, -38]))
    map.plot(xb, yb, 'k')

  # Save figure into output file...
  if (output):
    pp.savefig(output, bbox_inches='tight')

  # ...or show it
  else:
    pp.show()
    
#-------------------------------------------------------------------------------
# Routine : station_map
# Purpose : Draw map of stations
# Author  : P. Rebischung
# Created : 08-Feb-2012
#
# Changes :
#
# Input   : - lon         : Longitudes
#           - lat         : Latitudes
#           - code        : 4-char codes
#           - write_codes : Write 4-char codes names on map? Default is True.
#           - title       : Map title. Default is None.
#           - output      : Output file. Default is None.
# Output  : 
#-------------------------------------------------------------------------------
def station_map(lon, lat, code, write_codes=True, title=None, output=None):

  # Draw basemap
  pp.figure(figsize=(16, 12))
  map = Basemap(projection='robin', lon_0=0, resolution='c', area_thresh=50000.)
  map.drawmapboundary(fill_color=(0.82, 0.86, 1, 1)) 
  map.fillcontinents(color='white', lake_color=(0.82, 0.86, 1, 1), zorder=1)
  map.drawcoastlines(color='blue', linewidth=0.3)

  # Add title if necessary
  if (title):
    pp.title(title, position=(0.5, 1.06))

  # Plot points
  (x, y) = list(map(lon, lat))
  map.plot(x, y, '.', color='black', markersize=6)

  # Write 4-char codes if requested
  if (write_codes):
    for i in range(len(code)):
      pp.text(x[i], y[i], code[i], fontsize=14)

  # Save figure into output file...
  if (output):
    pp.savefig(output, bbox_inches='tight')

  # ...or show it
  else:
    pp.show()
    
#-------------------------------------------------------------------------------
# Routine : meanpole
# Purpose : Compute long term mean pole
# Author  : P. Rebischung
# Created : 20-Dec-2013
#
# Changes :
#
# Input   : - d      : MJDs
#           - meanpm : Keyword indicating which mean pole model should be used
#                      ('IERS2003' or 'IERS2010'). Default is 'IERS2010'.
# Output  : - xm     : Mean XPO (mas)
#           - ym     : Mean YPO (mas)
#-------------------------------------------------------------------------------
def meanpole(d, meanpm='IERS2010'):
  
  # Years / Years since 2000.0
  dy = (d - 51544.) / 365.25
  y  = dy + 2000.
  
  # IERS 2010 mean pole model
  if (meanpm == 'IERS2010'):
    
    # Initialize outputs
    xm = numpy.zeros(len(d))
    ym = numpy.zeros(len(d))
  
    # Mean pole before 2010.0
    ind = numpy.nonzero(y <= 2010.)[0]
    xm[ind] =  55.974 + 1.8243*dy[ind] + 0.18413*dy[ind]**2 + 0.007024*dy[ind]**3
    ym[ind] = 346.346 + 1.7896*dy[ind] - 0.10729*dy[ind]**2 - 0.000908*dy[ind]**3
    
    # Mean pole after 2010.0
    ind = numpy.nonzero(y > 2010.)[0]
    xm[ind] =  23.513 + 7.6141*dy[ind]
    ym[ind] = 358.891 - 0.6287*dy[ind]
    
  # IERS 2003 mean pole model
  elif (meanpm == 'IERS2003'):
    xm =  54. + 0.83*dy
    ym = 357. + 3.95*dy
  
  return (xm, ym)

#-------------------------------------------------------------------------------
# Routine : read_opoleload
# Purpose : Read grid of ocean pole tide loading coefficients
# Author  : P. Rebischung
# Created : 20-Dec-2013
#
# Changes :
#
# Input   : - inp   : Grid file
# Output  : - opole : Object containing interpolating splines for the 6 coefficients
#-------------------------------------------------------------------------------
def read_opoleload(inp):
  
  lon = numpy.arange(0.25, 360, 0.5)
  lat = numpy.arange(-89.75, 90, 0.5)
  (uhr, uhi, unr, uni, uer, uei) = numpy.loadtxt(inp, unpack=True, skiprows=14, usecols=list(range(2, 8)))

  opole = record()

  opole.uhr = interp.RectBivariateSpline(lon, lat, uhr.reshape((360, 720)).T)
  opole.uhi = interp.RectBivariateSpline(lon, lat, uhi.reshape((360, 720)).T)
  opole.unr = interp.RectBivariateSpline(lon, lat, unr.reshape((360, 720)).T)
  opole.uni = interp.RectBivariateSpline(lon, lat, uni.reshape((360, 720)).T)
  opole.uer = interp.RectBivariateSpline(lon, lat, uer.reshape((360, 720)).T)
  opole.uei = interp.RectBivariateSpline(lon, lat, uei.reshape((360, 720)).T)
  
  return opole

#-------------------------------------------------------------------------------
# Routine : compute_opoleload
# Purpose : Compute ocean pole tide loading displacements
# Author  : P. Rebischung
# Created : 20-Dec-2013
#
# Changes :
#
# Input   : - opole  : Object containing interpolating splines for the 6 coefficients
#           - lon    : Station longitudes (size n)
#           - lat    : Station latitudes  (size n)
#           - mjd    : MJDs (size t)
#           - xpo    : X pole coordinates
#           - ypo    : Y pole coordinates
#           - meanpm : Keyword indicating which mean pole model should be used
#                      ('IERS2003' or 'IERS2010'). Default is 'IERS2010'.
# Output  : - denh   : Station displacements (size [n] x [t] x 3)
#-------------------------------------------------------------------------------
def compute_opoleload(opole, lon, lat, mjd, xpo, ypo, meanpm='IERS2010'):
  
  # If only one station,
  if (type(lon).__name__ == 'float'):
    lon = numpy.array([lon])
    lat = numpy.array([lat])
    
  # If only one date,
  if (type(mjd).__name__ == 'float'):
    mjd = numpy.array([mjd])

  # Get number of stations - Add 360 degrees to negative longitudes
  n = len(lon)
  ind = numpy.nonzero(lon < 0)[0]
  lon[ind] = lon[ind] + 360

  # Get number of epochs
  t = len(mjd)
  
  # Scaling constants
  Hp = sqrt(8*pi/15) * Oe**2 * ae**4 / GM
  K  = 4*pi*G*ae*1025*Hp / (3*ge)
  
  # Evaluate ocean pole tide loading coefficients at requested stations
  uhr = opole.uhr.ev(lon, lat)
  uhi = opole.uhi.ev(lon, lat)
  unr = opole.unr.ev(lon, lat)
  uni = opole.uni.ev(lon, lat)
  uer = opole.uer.ev(lon, lat)
  uei = opole.uei.ev(lon, lat)

  # Compute mean pole
  (xpm, ypm) = meanpole(mjd, meanpm)
  
  # Compute wobble parameters
  m1 =  (xpo - xpm) * mas2rad
  m2 = -(ypo - ypm) * mas2rad
 
  # Compute displacement amplitudes
  ar = K * (0.6870*m1 + 0.0036*m2)
  ai = K * (0.6870*m2 - 0.0036*m1)
    
  # Compute displacements
  denh = numpy.zeros((n, t, 3))
  for i in range(n):
    denh[i,:,0] = ar*uer[i] + ai*uei[i]
    denh[i,:,1] = ar*unr[i] + ai*uni[i]
    denh[i,:,2] = ar*uhr[i] + ai*uhi[i]
      
  return numpy.squeeze(denh)

#-------------------------------------------------------------------------------
# Routine : compute_psd
# Purpose : Compute post-seismic deformations and associated uncertainties
# Author  : P. Rebischung
# Created : 21-Aug-2015
#
# Changes :
#
# Input   : - snx : sinex object containing post-seismic deformation models
#           - sta : 4-char station ID
#           - t   : date in SINEX format
# Output  : - dx  : E,N,H post-seismic deformations at time t [m]
#           - sx  : associated formal errors [m]
#-------------------------------------------------------------------------------
def compute_psd(snx, sta, t):
  
  # Initializations
  dx = numpy.zeros(3)
  sx = numpy.zeros(3)
  mjd = date.from_snxepoch(t).mjd
   
  # Get indices of parameters describing the East component of the post-seismic deformations
  ind = []
  for i in range(snx.npar):
    p = snx.param[i]
    if ((p.code == sta) and (p.type[5] == 'E') and (earlier(p.tref, t))):
      ind.append(i)
      
  # Loop over model functions
  A = numpy.zeros(len(ind))
  for i in range(0, len(ind), 2):
    
    # Useful things
    mjd0 = date.from_snxepoch(snx.param[ind[i]].tref).mjd
    dt = (mjd - mjd0) / 365.25
    amp = snx.x[ind[i]]
    tau = snx.x[ind[i+1]]
    
    # Case of an exponential
    if (snx.param[ind[i]].type[1:4] == 'EXP'):
      
      # Compute associated deformation
      da = 1. - exp(-dt / tau)
      dx[0] = dx[0] + amp * da
      
      # Compute partial derivatives of the model parameters
      A[i] = da
      A[i+1] = -amp * dt * (1-da) / tau**2
      
    # Case of a logarithm
    elif (snx.param[ind[i]].type[1:4] == 'LOG'):
      
      # Compute associated deformation
      da = log(1 + dt / tau)
      dx[0] = dx[0] + amp * da
      
      # Compute partial derivatives of the model parameters
      A[i] = da
      A[i+1] = -amp * dt / (1 + dt / tau) / tau**2
      
  # Compute formal error of the East component of the post-seismic deformations
  if (len(ind) > 0):
    Q = snx.Q[numpy.ix_(ind,ind)]
    sx[0] = sqrt(dot(A, dot(Q, A.T)))

  # Get indices of parameters describing the North component of the post-seismic deformations
  ind = []
  for i in range(snx.npar):
    p = snx.param[i]
    if ((p.code == sta) and (p.type[5] == 'N') and (earlier(p.tref, t))):
      ind.append(i)
      
  # Loop over model functions
  A = numpy.zeros(len(ind))
  for i in range(0, len(ind), 2):
    
    # Useful things
    mjd0 = date.from_snxepoch(snx.param[ind[i]].tref).mjd
    dt = (mjd - mjd0) / 365.25
    amp = snx.x[ind[i]]
    tau = snx.x[ind[i+1]]
    
    # Case of an exponential
    if (snx.param[ind[i]].type[1:4] == 'EXP'):
      
      # Compute associated deformation
      da = 1. - exp(-dt / tau)
      dx[1] = dx[1] + amp * da
      
      # Compute partial derivatives of the model parameters
      A[i] = da
      A[i+1] = -amp * dt * (1-da) / tau**2
      
    # Case of a logarithm
    elif (snx.param[ind[i]].type[1:4] == 'LOG'):
      
      # Compute associated deformation
      da = log(1 + dt / tau)
      dx[1] = dx[1] + amp * da
      
      # Compute partial derivatives of the model parameters
      A[i] = da
      A[i+1] = -amp * dt / (1 + dt / tau) / tau**2
      
  # Compute formal error of the North component of the post-seismic deformations
  if (len(ind) > 0):
    Q = snx.Q[numpy.ix_(ind,ind)]
    sx[1] = sqrt(dot(A, dot(Q, A.T)))

  # Get indices of parameters describing the Up component of the post-seismic deformations
  ind = []
  for i in range(snx.npar):
    p = snx.param[i]
    if ((p.code == sta) and (p.type[5] == 'H') and (earlier(p.tref, t))):
      ind.append(i)
      
  # Loop over model functions
  A = numpy.zeros(len(ind))
  for i in range(0, len(ind), 2):
    
    # Useful things
    mjd0 = date.from_snxepoch(snx.param[ind[i]].tref).mjd
    dt = (mjd - mjd0) / 365.25
    amp = snx.x[ind[i]]
    tau = snx.x[ind[i+1]]
    
    # Case of an exponential
    if (snx.param[ind[i]].type[1:4] == 'EXP'):
      
      # Compute associated deformation
      da = 1. - exp(-dt / tau)
      dx[2] = dx[2] + amp * da
      
      # Compute partial derivatives of the model parameters
      A[i] = da
      A[i+1] = -amp * dt * (1-da) / tau**2
      
    # Case of a logarithm
    elif (snx.param[ind[i]].type[1:4] == 'LOG'):
      
      # Compute associated deformation
      da = log(1 + dt / tau)
      dx[2] = dx[2] + amp * da
      
      # Compute partial derivatives of the model parameters
      A[i] = da
      A[i+1] = -amp * dt / (1 + dt / tau) / tau**2
      
  # Compute formal error of the Up component of the post-seismic deformations
  if (len(ind) > 0):
    Q = snx.Q[numpy.ix_(ind,ind)]
    sx[2] = sqrt(dot(A, dot(Q, A.T)))

  return (dx, sx)