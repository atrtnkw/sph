from yt.mods import *
import sys

filetype = "cdet_hdf5_chk_"

argvs = sys.argv
argc  = len(argvs)
if (argc != 5):
    print ('python plot1d.py <idir> <tbgn> <tend> <dtsp>')
    quit()

idir = argvs[1]
tbgn = int(argvs[2])
tend = int(argvs[3])
dtsp = int(argvs[4])

for time in range(tbgn,tend,dtsp):
    ctime = str(time).zfill(4)
    ifile = idir + filetype + ctime
    ofile = "mesh_t" + ctime + ".dat"
    
    if not os.path.exists(ifile) or os.path.exists(ofile) :
        continue

    ds = load(ifile)
    ad = ds.h.all_data()
    obj = open(ofile, 'w')
    for i in range(ad["temp"].size):
        print (ad["r"][i], ad["dens"][i], ad["temp"][i], ad["eint"][i], ad["pres"][i], ad["velx"][i], ad["enuc"][i], ad["he4 "][i], ad["c12 "][i], ad["o16 "][i], ad["ne20"][i], ad["mg24"][i], ad["si28"][i], ad["s32 "][i], ad["ar36"][i], ad["ca40"][i], ad["ti44"][i], ad["cr48"][i], ad["fe52"][i], ad["ni56"][i], file=obj)
