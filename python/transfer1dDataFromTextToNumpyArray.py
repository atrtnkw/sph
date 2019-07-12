import sys
import numpy as np

argv = sys.argv
argc = len(argv)

if (argc != 4):
    print ('python xxx.py <ifile> <icolumn> <otype>')
    quit()

ifile = argv[1]
icolm = int(argv[2])
otype = argv[3]
ofile = argv[3] + ".npy"

dens = []
#for line in open('input.data').readlines():
for line in open(ifile).readlines():
    data = line.split()
    dens.append(float(data[icolm]))

dens = np.array(dens,dtype='float64')
#dens = np.array(dens,dtype='float32')

np.save(otype, dens)
print(np.load(ofile))
