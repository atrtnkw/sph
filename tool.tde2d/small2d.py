from yt.mods import *
import sys

filetype = "tde2d_hdf5_chk_"
argvs = sys.argv
argc  = len(argvs)
if (argc != 4):
    print ('python plot1d.py <tbgn> <tend> <dtsp>')
    quit()

tbgn = int(argvs[1])
tend = int(argvs[2])
dtsp = int(argvs[3])

for time in range(tbgn,tend,dtsp):
    ctime = str(time).zfill(4)
    ifile = filetype + ctime
    ofile = "small_t" + ctime + ".dat"

    if (not os.path.exists(ifile)) or os.path.exists(ofile) :
        continue

    ds = load(ifile)
    ad = ds.h.all_data()
    obj = open(ofile, 'w')
    length = range(ad["temp"].size)
    for i in length:
        print (ad["x"][i], ad["y"][i], ad["dens"][i], ad["temp"][i], file=obj)
