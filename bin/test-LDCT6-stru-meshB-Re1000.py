import matplotlib.pyplot as plt
import numpy as np
import sys
from pylab import *
import os


#inpfile   =  sys.argv[1]

inpfile = "LDCT6-stru-meshB"
conv_data_ref = "convergence-data-LDCT6-stru-meshB-Re1000.dat"

print("Started the simulation \n")
#print("Check 'simulation.log' for progress \n")

#cmd="./incexplicitSerial  " + inpfile + "  >  simulation.log"
cmd="./incexplicitSerial  " + inpfile

os.system(cmd)

print("Completed the simulation \n")


xR=np.loadtxt(conv_data_ref)
x1=np.loadtxt("convergence-data.dat")

xDiff = xR[:,2] - x1[:,2]

errnorm = np.linalg.norm(xDiff)

print("\nError norm = %12.6f \n" % errnorm)

if( errnorm < 1.0e-10):
    print("Test for '%s' is SUCCESSFUL \n" % inpfile)
else:
    print("Test for '%s' is NOT SUCCESSFUL \n" % inpfile)

plt.plot(xR[:,1],np.log10(xR[:,2]), 'r', label="Reference", linewidth=2.0, markersize=8.0)
plt.plot(x1[:,1],np.log10(x1[:,2]), 'k.', label="Present", linewidth=2.0, markersize=8.0)
plt.legend(loc="upper right")
plt.grid('on')


plt.savefig("convergence.png", dpi=500)

