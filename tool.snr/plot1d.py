import sys
from yt.mods import *

filetype = "sedov_hdf5_chk_"

argvs = sys.argv
argc  = len(argvs)
if (argc != 4):
#    print 'python plot1d.py <tbgn> <tend> <dtsp>'
    print ('python plot1d.py <tbgn> <tend> <dtsp>')
    quit()

tbgn = int(argvs[1])
tend = int(argvs[2])
dtsp = int(argvs[3])

for time in range(tbgn,tend,dtsp):
    ctime = str(time).zfill(4)
    ifile = filetype + ctime
    ofile = "mesh_t" + ctime + ".dat"

    if not os.path.exists(ifile) or os.path.exists(ofile) :
        continue

    ds = load(ifile)
    ad = ds.h.all_data()
    obj = open(ofile, 'w')
    for i in range(ad["temp"].size):
        print ("%e %e %e %e %e %e" % (ad["r"][i], ad["dens"][i], ad["temp"][i], ad["pres"][i], ad["ener"][i], ad["velx"][i]), file=obj)
#        print ("%e %e %e %e %e %e" % (ad["x"][i], ad["dens"][i], ad["temp"][i], ad["pres"][i], ad["ener"][i], ad["velx"][i]), file=obj)

