#!/usr/bin/env python

#### Script to calculate statitistics of SMS scan

import os,sys,math
import numpy as np
import matplotlib.pyplot as plt

# Number of events: min(xfactor*xsec*ifb, maxEvents) (always in thousands)
ifb, xfactor, maxEvents, minLumi = 20, 40, 1000, 40
small_ystep = 50
font_size = 8

model = "T5qqqqVV"
model = "T1bbbb"
model = "T1tttt"
model = "T2tt"

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
minEvents = 10
diagStep = 50
maxDM = 1000
addDiag = []
if(model=="T2tt"):
  font_size = 7
  diagStep = 25
  maxDM = 250
  minLumi = 1e-20
  addDiag = [183, 167, 85]
  scanBlocks.append(gridBlock(100,  701, 50, 50, maxDM, diagStep, minEvents))
  minDM, ymed, ymax = 75, 0, 600 
  
if(model=="T1tttt"):
  grid_type = "2016"
  if(grid_type == "2015"):
    scanBlocks.append(gridBlock(600,  1000, 100, 100, 425, 25, 50))
    scanBlocks.append(gridBlock(1000, 1600, 50, 50, 425, 25, 100))
    scanBlocks.append(gridBlock(1600, 2001, 50, 50, 425, 25, 50))
    minDM, ymed, ymax = 225, 300, 1450 
  else:
    minEvents = 20
    scanBlocks.append(gridBlock(600,  1200, 100, 100, maxDM, diagStep, minEvents))
    scanBlocks.append(gridBlock(1200, 2000, 50, 100, maxDM, diagStep, minEvents))
    scanBlocks.append(gridBlock(2000, 2301, 50, 100, maxDM, diagStep, minEvents))
    minDM, ymed, ymax = 225, 300, 1600 

if(model=="T1bbbb"):
  scanBlocks.append(gridBlock(600,  1200, 100, 100, maxDM, diagStep, minEvents))
  scanBlocks.append(gridBlock(1200, 2000, 50, 100, maxDM, diagStep, minEvents))
  scanBlocks.append(gridBlock(2000, 2301, 50, 100, maxDM, diagStep, minEvents))
  minDM, ymed, ymax = 25, 500, 1600 

if(model=="T5qqqqVV"):
  scanBlocks.append(gridBlock(600,  1200, 100, 100, maxDM, diagStep, minEvents))
  scanBlocks.append(gridBlock(1200, 2000, 50, 100, maxDM, diagStep, minEvents))
  scanBlocks.append(gridBlock(2000, 2301, 50, 100, maxDM, diagStep, minEvents))
  minDM, ymed, ymax = 125, 400, 1600 

# Fit to gluino cross-section in fb
def xsecGlu(mass):
  return 4.563e+17*math.pow(mass, -4.761*math.exp(5.848e-05*mass))
  #return 2.02584e+17*math.pow(mass, -4.6368*math.exp(6.16673e-05*mass)) # Fit to [1000, 2000]

# Fit to stop cross-section in fb
def xsecStop(mass):
  if mass < 300: return 319925471928717.38*math.pow(mass, -4.10396285974583*math.exp(mass*0.0001317804474363))
  else: return 6953884830281245*math.pow(mass, -4.7171617288678069*math.exp(mass*6.1752771466190749e-05))

# Number of events for mass point, in thousands
def events(mass):
  if("T2" in model): xsec = xsecStop(mass)
  else: xsec = xsecGlu(mass)
  nev = min(xfactor*xsec*ifb, maxEvents*1000)
  if nev < xsec*minLumi: nev = xsec*minLumi
  nev = max(nev/1000, minEvents)
  return math.ceil(nev) # Rounds up

def makePlot(mpoints, type):
  plt.figure(figsize=(17,10))
  if("T1" in model or "T5" in model): plt.xlabel('$m(\widetilde{g})$ [GeV]', fontsize=18)
  if("T2tt"==model): plt.xlabel('$m(\widetilde{t})$ [GeV]', fontsize=18)
  plt.ylabel('$m(\chi^0_1)$ [GeV]', fontsize=18)
  Ntot = 0
  for col in mpoints:
    for mpoint in col:
      if("T2" in model): xsec = xsecStop(mpoint[0])
      else: xsec = xsecGlu(mpoint[0])
      nev = mpoint[2]
      Ntot += nev
      if type == 'events': val = nev
      if type == 'factor': val = nev/xsec/ifb*1000
      if type == 'lumi': val = nev/xsec*1000
      if val<1000:
        plt.text(mpoint[0],mpoint[1], "{0:.0f}".format(val), 
                 verticalalignment='center', horizontalalignment='center', fontsize=font_size)
      else:
        plt.text(mpoint[0],mpoint[1], "{0:.1f}".format(val/1000), fontweight='bold', 
                 verticalalignment='center', horizontalalignment='center', fontsize=font_size, color='red')
  plt.axis([xmin-10, xmax+10, -50, ymax+50])
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
    begDiag = max(ymed, mx-block.maxDM)
    # Adding bulk points
    if (mx-block.xmin)%block.xstep == 0 :
      for my in range(ymin, begDiag, block.ystep):
        if my > ymax: continue
        for dm in addDiag:
          dm_before = mx
          if(len(col)>0): dm_before = col[-1][0]-col[-1][1]
          dm_after = mx - my
          if(dm<dm_before and dm>dm_after):
            nev = events(mx)
            col.append([mx,mx-dm, nev])
            Nbulk += nev
            
        nev = events(mx)
        col.append([mx,my, nev])
        Nbulk += nev
    # Adding diagonal points
    for my in range(begDiag, mx-minDM+1, block.dstep):
      if my > ymax: continue
      for dm in addDiag:
        dm_before = mx
        if(len(col)>0): dm_before = col[-1][0]-col[-1][1]
        dm_after = mx - my
        #print "dm "+str(dm)+", bef "+str(dm_before)+", aft "+str(dm_after)
        if(dm<dm_before and dm>dm_after):
          nev = events(mx)
          col.append([mx,mx-dm, nev])
          Ndiag += nev

      nev = events(mx)
      col.append([mx,my, nev])
      Ndiag += nev
    if(my !=  mx-minDM and mx-minDM <= ymax):
      my = mx-minDM
      for dm in addDiag:
        dm_before = mx
        if(len(col)>0): dm_before = col[-1][0]-col[-1][1]
        dm_after = mx - my
        print "dm "+str(dm)+", bef "+str(dm_before)+", aft "+str(dm_after)
        if(dm<dm_before and dm>dm_after):
          nev = events(mx)
          col.append([mx,mx-dm, nev])
          Ndiag += nev

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
