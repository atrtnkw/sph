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
    ofile = "mesh_t" + ctime + ".dat" + ".xy"

    if not os.path.exists(ifile) or os.path.exists(ofile) :
        continue

    ds = load(ifile)
#    ad = ds.r[:,:,0.00390625]
    ad = ds.r[:,:,3.90625e16]
    obj = open(ofile, 'w')
    for i in range(ad["temp"].size):
        print ("%e %e %e %e %e %e %e" % (ad["x"][i], ad["y"][i], ad["z"][i], ad["dens"][i], ad["temp"][i], ad["pres"][i], ad["ener"][i]), file=obj)
