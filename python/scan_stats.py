#!/usr/bin/env python

#### Script to calculate statitistics of SMS scan

import os,sys,math
import matplotlib.pyplot as plt

# Number of events: min(xfactor*xsec*ifb, maxEvents) (always in thousands)
ifb, xfactor, maxEvents = 20, 20, 150

# Parameters to define the bulk of the grid
xaxis = [[800, 1200, 100], [1200, 2201, 50]] #[[xmin, xmax, xstep]]
ymin, ymax, ystep = 0, 1500, 100 
# Parameters for grid near diagonal
minDM, maxDM = 225, 400
dstep = 25 # needs to be smaller than all xsets in xaxis


# Fit to gluino cross-section in fb
def xsecGlu(mass):
  return 2.02584e+17*math.pow(mass, -4.6368*math.exp(6.16673e-05*mass))

# Number of events for mass point, in thousands
def events(mass):
  nev = min(xfactor*xsecGlu(mass)*ifb, maxEvents*1000)
  return math.ceil(nev/1000) # Rounds up

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
      plt.text(mpoint[0],mpoint[1], "{0:.0f}".format(val), 
               verticalalignment='center', horizontalalignment='center', fontsize=7)
  plt.axis([xmin-100, xmax+100, -50, ymax+100])
  plt.grid(True)
  if type == 'events': title = 'Thousands of events to generate'
  if type == 'factor': title = 'Times more MC events than expected for '+str(ifb)+' fb$^{-1}$'
  tot_s = ' ('+"{0:.1f}".format(Ntot/1000)+' million events in the scan)'
  plt.title(title+tot_s, fontweight='bold')
  pname = 'sms_'+type+'.pdf'
  plt.savefig(pname)
  print ' open '+pname
  return Ntot


mpoints = []
xmin, xmax = 9999, 0
for xpar in xaxis:
  for mx in range(xpar[0], xpar[1], dstep):
    xmin = min(xmin, xpar[0])
    xmax = max(xmax, xpar[1])
    col = []
    if (mx-xpar[0])%xpar[2] == 0 :
      for my in range(ymin, mx-maxDM, ystep):
        if my > ymax: continue
        col.append([mx,my, events(mx)])
    for my in range(mx-maxDM, mx-minDM+1, dstep):
      if my > ymax: continue
      col.append([mx,my, events(mx)])
    mpoints.append(col)


makePlot(mpoints, 'events')
Ntot = makePlot(mpoints, 'factor')

print '\nTotal scan amounts to '+"{0:.1f}".format(Ntot/1000)+" million events\n"
