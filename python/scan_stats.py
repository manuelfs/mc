#!/usr/bin/env python

#### Script to calculate statitistics of SMS scan

import os,sys,math
import matplotlib.pyplot as plt

# Number of events: min(xfactor*xsec*ifb, maxEvents) (always in thousands)
ifb, xfactor, maxEvents = 30, 20, 500

# Parameters to define the bulk of the grid
xaxis = [[800, 1200, 100], [1200, 2201, 50]] #[[xmin, xmax, xstep]]
ymin, ystep = 0, 100 
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

mpoints = []
Ntot = 0
xmin, xmax = 9999, 0
for xpar in xaxis:
  for mx in range(xpar[0], xpar[1], dstep):
    xmin = min(xmin, xpar[0])
    xmax = max(xmax, xpar[1])
    col = []
    if (mx-xpar[0])%xpar[2] == 0 :
      for my in range(ymin, mx-maxDM, ystep):
        nev = events(mx)
        col.append([mx,my, nev])
        Ntot += nev
    for my in range(mx-maxDM, mx-minDM+1, dstep):
      nev = events(mx)
      col.append([mx,my, nev])
      Ntot += nev
    mpoints.append(col)

print str(Ntot/1000)+" million events"


plt.xlabel('$m(\widetilde{g})$ [GeV]')
plt.ylabel('$m(\chi^0_1)$ [GeV]')
plt.title('Total number of events in scan: '+str(Ntot/1000)+" million")
for col in mpoints:
  for mpoint in col:
    plt.text(mpoint[0],mpoint[1], "{0:.0f}".format(mpoint[2]), 
             verticalalignment='center', horizontalalignment='center', fontsize=10)
plt.axis([xmin-100, xmax+100, -50, xmax-minDM+100])
plt.grid(True)
plt.show()

