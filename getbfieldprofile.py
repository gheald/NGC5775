#!/usr/bin/env python

from numpy import *
from pylab import *
import pyfits
import argparse

def main(args):

	R = pyfits.getdata('galR.fits')
	z = pyfits.getdata('galz.fits')
	B = pyfits.getdata('n5775_bfield_equipart.fits')

	# Here are the quadrant definitions.
	# Quadr	Sign(R)	Sign(z)
	# NE	-	-
	# NW	-	+
	# SE	+	-
	# SW	+	+

	q = args.quadrant
	s = { 'NE': (-1.,-1.),
	      'NW': (-1.,1.),
	      'SE': (1.,-1.),
	      'SW': (1.,1.)   }

	"""
	THIS IS THE OLD STUFF

	#gp = logical_and(R>-0.035,R<0.035)
	#gp = logical_and(gp,z>-0.00375)
	#gp = logical_and(gp,z<0.00125)

	#plot(R[gp]*3600./206265.*24800.,B[gp],'ko')
	#xlabel('R (kpc)')
	#ylabel('B (uG)')
	#text(10.,25.,'|z| < %.2f kpc'%(0.0025*3600./206265.*24800.),size='large')
	#savefig('n5775_midplaneBeq.png',dpi=200,bbox_inches='tight')
	#show()

	"""

	ny = B.shape[0]
        nx = B.shape[1]
	dc = pi/180.*24800. # distance conversion (degr to kpc)
	zaxis = arange(0.,11.1,1.)/dc
	Rzgrid = zeros((1,len(zaxis)))
	ccgrid = zeros((1,len(zaxis)))
	ergrid = zeros((1,len(zaxis)))
	for x in range(nx):
        	for y in range(ny):
			# select quadrant
			if s[q][0]*R[y,x] < 0. or s[q][1]*z[y,x] < 0.: continue
			# reject pixels outside of region
                	if abs(R[y,x])>90./3600. or abs(z[y,x])>11./dc: continue
                	#if R[y,x]>0./3600. or z[y,x]<0./dc: continue
                	if isnan(B[y,x]): continue
                	ri = 0#where(abs(Raxis)-abs(R[y,x])==min(abs(Raxis)-abs(R[y,x])))
                	#zi = where(abs(zaxis-abs(z[y,x]))==min(abs(zaxis-abs(z[y,x]))))
                	zi = where(abs(zaxis-z[y,x])==min(abs(zaxis-z[y,x])))
                	Rzgrid[ri,zi] += B[y,x]
                	ccgrid[ri,zi] += 1
	Rzgrid[ccgrid!=0] /= ccgrid[ccgrid!=0]
	for x in range(nx):
        	for y in range(ny):
			# select quadrant
			if s[q][0]*R[y,x] < 0. or s[q][1]*z[y,x] < 0.: continue
			# reject pixels outside of region
                	if abs(R[y,x])>90./3600. or abs(z[y,x])>11./dc: continue
                	#if R[y,x]>0./3600. or z[y,x]<0./dc: continue
                	if isnan(B[y,x]): continue
                	ri = 0#where(abs(Raxis)-abs(R[y,x])==min(abs(Raxis)-abs(R[y,x])))
                	#zi = where(abs(zaxis-abs(z[y,x]))==min(abs(zaxis-abs(z[y,x]))))
                	zi = where(abs(zaxis-z[y,x])==min(abs(zaxis-z[y,x])))
                	ergrid[ri,zi] += (B[y,x]-Rzgrid[ri,zi])**2
	ergrid[ccgrid!=0] /= ccgrid[ccgrid!=0]
	ergrid = sqrt(ergrid)
	outputfile = open(args.outfile,'w')
	for i in range(len(zaxis)):
        	print >>outputfile, abs(zaxis[i]*dc*1000.), Rzgrid[0,i], ergrid[0,i] #args.noise/sqrt(ccgrid[0,i])
	outputfile.close()

ap = argparse.ArgumentParser()
ap.add_argument('--outfile','-o',help='Output file [default Bprof.txt]',default='Bprof.txt')
ap.add_argument('quadrant',help='Specify quadrant (choose from NE,NW,SE,SW)',choices=('NE','NW','SE','SW'))
args = ap.parse_args()
main(args)

