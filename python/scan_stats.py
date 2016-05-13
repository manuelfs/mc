#!/usr/bin/env python

#### Script to calculate statitistics of SMS scan

import os,sys,math
import numpy as np
import matplotlib.pyplot as plt

# Number of events: min(xfactor*xsec*ifb, maxEvents) (always in thousands)
ifb, xfactor, maxEvents, minLumi = 20, 40, 150, 40
small_ystep = 50

model = "T5qqqqVV"
model = "T1bbbb"
model = "T1tttt"

# Parameters that define the grid in the bulk and diagonal
class gridBlock:
  def __init__(self, xmin, xmax, xstep, ystep, maxDM, dstep, minEvents):
    self.xmin = xmin
    self.xmax = xmax
    self.xstep = xstep
    self.ystep = ystep
    self.maxDM = maxDM
    self.dstep = dstep
    self.minEvents = minEvents

# Parameters to define the grid of the various models
ymin = 0
scanBlocks = []
grid_type = "2016"
if(model=="T1tttt"):
  if(grid_type == "2015"):
    scanBlocks.append(gridBlock(600,  1000, 100, 100, 425, 25, 50))
    scanBlocks.append(gridBlock(1000, 1600, 50, 50, 425, 25, 100))
    scanBlocks.append(gridBlock(1600, 2001, 50, 50, 425, 25, 50))
    minDM, ymed, ymax = 225, 300, 1450 
  else:
    scanBlocks.append(gridBlock(600,  1200, 100, 100, 1000, 50, 10))
    scanBlocks.append(gridBlock(1200, 2000, 50, 100, 1000, 50, 10))
    scanBlocks.append(gridBlock(2000, 2301, 50, 100, 1000, 50, 5))
    minDM, ymed, ymax = 225, 300, 1600 

if(model=="T1bbbb"):
  scanBlocks.append(gridBlock(600,  1200, 100, 100, 1000, 50, 10))
  scanBlocks.append(gridBlock(1200, 2000, 50, 100, 1000, 50, 10))
  scanBlocks.append(gridBlock(2000, 2301, 50, 100, 1000, 50, 5))
  minDM, ymed, ymax = 25, 500, 1600 

if(model=="T5qqqqVV"):
  scanBlocks.append(gridBlock(600,  1200, 100, 100, 1000, 50, 10))
  scanBlocks.append(gridBlock(1200, 2000, 50, 100, 1000, 50, 10))
  scanBlocks.append(gridBlock(2000, 2301, 50, 100, 1000, 50, 5))
  minDM, ymed, ymax = 125, 400, 1600 

# Fit to gluino cross-section in fb
def xsecGlu(mass):
  return 2.02584e+17*math.pow(mass, -4.6368*math.exp(6.16673e-05*mass))

# Number of events for mass point, in thousands
def events(mass):
  xsec = xsecGlu(mass)
  nev = min(xfactor*xsec*ifb, maxEvents*1000)
  if nev < xsec*minLumi: nev = xsec*minLumi
  nev = max(nev/1000, minEvents)
  return math.ceil(nev) # Rounds up

def makePlot(mpoints, type):
  plt.figure(figsize=(17,10))
  plt.xlabel('$m(\widetilde{g})$ [GeV]', fontsize=18)
  plt.ylabel('$m(\chi^0_1)$ [GeV]', fontsize=18)
  Ntot = 0
  for col in mpoints:
    for mpoint in col:
      nev = mpoint[2]
      Ntot += nev
      if type == 'events': val = nev
      if type == 'factor': val = nev/xsecGlu(mpoint[0])/ifb*1000
      if type == 'lumi': val = nev/xsecGlu(mpoint[0])*1000
      plt.text(mpoint[0],mpoint[1], "{0:.0f}".format(val), 
               verticalalignment='center', horizontalalignment='center', fontsize=8)
  plt.axis([xmin-100, xmax+100, -50, ymax+100])
  plt.xticks(np.arange(xmin, xmax, 200))
  plt.yticks(np.arange(ymin, ymax+100, 200))
  plt.grid(True)
  if type == 'events': title = 'Thousands of '+model+' events to generate'
  if type == 'factor': title = 'Times more '+model+' MC events than expected in data for '+str(ifb)+' fb$^{-1}$'
  if type == 'lumi': title = 'Equivalent '+model+' MC luminosity in fb$^{-1}$'
  tot_s = ' ('+"{0:.1f}".format(Ntot/1000)+' million events in the scan)'
  plt.title(title+tot_s, fontweight='bold')
  pname = model+'_'+type+'.pdf'
  plt.savefig(pname, bbox_inches='tight')
  print ' open '+pname
  return Ntot


mpoints = []
Nevents = []
xmin, xmax = 9999, 0
for block in scanBlocks:
  Nbulk, Ndiag = 0, 0
  minEvents = block.minEvents
  for mx in range(block.xmin, block.xmax, block.dstep):
    xmin = min(xmin, block.xmin)
    xmax = max(xmax, block.xmax)
    col = []
    my = 0
    begDiag = max(ymed, mx-block.maxDM,)
    # Adding bulk points
    if (mx-block.xmin)%block.xstep == 0 :
      for my in range(ymin, begDiag, block.ystep):
        if my > ymax: continue
        nev = events(mx)
        col.append([mx,my, nev])
        Nbulk += nev
    # Adding diagonal points
    for my in range(begDiag, mx-minDM+1, block.dstep):
      if my > ymax: continue
      nev = events(mx)
      col.append([mx,my, nev])
      Ndiag += nev
    if(my !=  mx-minDM and mx-minDM <= ymax):
      my = mx-minDM
      nev = events(mx)
      col.append([mx,my, nev])
      Ndiag += nev
    mpoints.append(col)
  Nevents.append([Nbulk, Ndiag])


makePlot(mpoints, 'events')
Ntot = makePlot(mpoints, 'lumi')
#makePlot(mpoints, 'factor')

Ntot = Ntot/1e3
print '\nTotal scan amounts to '+"{0:.1f}".format(Ntot)+" million events\n"
for ind in range(len(scanBlocks)):
  Nbulk, Ndiag = Nevents[ind][0]/1e3, Nevents[ind][1]/1e3
  Nblock = Nbulk+Ndiag
  print "From "+'{:>4}'.format(scanBlocks[ind].xmin)+" to "+str(scanBlocks[ind].xmax)+": ",
  print "{0:>4.1f}".format(Nblock)+"M ("+"{0:>4.1f}".format(Nblock/Ntot*100)+" %) events, "+"{0:>4.1f}".format(Nbulk),
  print "M ("+"{0:>4.1f}".format(Nbulk/Ntot*100)+" %) in the bulk, "+"{0:>4.1f}".format(Ndiag)+"M (",
  print "{0:.1f}".format(Ndiag/Ntot*100)+" %) in the diagonal"

print
