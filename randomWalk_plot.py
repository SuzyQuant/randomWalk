#!/usr/bin/env python

import matplotlib
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

data = np.loadtxt('positionVSMeanSquareOutput.dat')
time = data[:, 0] + 1
MSDth = time
MSDnum = data[:, 1]
MSDerr = data[:, 2]
Fsq1 = data[:, 3]
Fsq1err = data[:, 6]


with PdfPages('positionVSMeanSquareOutput.pdf') as pdf:


	plt.xlim([1,500])
#	plt.ylim([,1.9])
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$MSD(t)$', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')
	
	plt.plot( time , MSDnum, 'r-*',label= r'Numeric', linewidth=1)
	plt.plot( time , MSDth, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======================	
	plt.xlim([1,500])
#	plt.ylim([,1.9])
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$MSD(t)$', fontsize=10)
	plt.xscale('log')
	plt.yscale('log')
	
	plt.errorbar(time , MSDnum, MSDerr, color='r', label= r'Numeric', linewidth=1)
	plt.plot( time , MSDth, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()

# ======================	
	plt.xlim([1,500])
#	plt.ylim([,1.9])
	plt.xlabel(r'$t$', fontsize=15)
	plt.ylabel(r'$MSD(t)$', fontsize=10)
	plt.xscale('linear')
	plt.yscale('linear')
	
	plt.errorbar(time , MSDnum, MSDerr)
#	plt.errorbar(time , MSDnum, MSDerr*10)	multiplied by 10 to increase error bar
	plt.plot( time , MSDth, 'c-',label= r'Theory', linewidth=1)

	plt.legend(loc=4)
	#plt.show()
	pdf.savefig()
	plt.close()