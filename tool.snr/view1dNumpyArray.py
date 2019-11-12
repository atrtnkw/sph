import sys
import numpy as np

argv = sys.argv
argc = len(argv)

if (argc != 2):
    print ('python xxx.py <ifile>')
    quit()

ifile = argv[1]

a = np.load(ifile)
np.set_printoptions(threshold=np.inf)
#np.set_printoptions(threshold=a.size)
#np.set_printoptions(threshold=128)
print(a)
