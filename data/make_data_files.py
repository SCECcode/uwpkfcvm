#!/usr/bin/env python

##
#  Builds the data files in the expected format from 
#
# from >> XX??  YY?? depth(KM) Vp/Vs (km/s)
#         -65.0 -70.0  0.0  3.17
#
# create binary vp.dat and vs.dat
#
## create a combined vs/vp and no density file(to be calculated)

import getopt
import sys
import subprocess
import struct
import array
import ssl
import certifi

if sys.version_info.major >= (3) :
  from urllib.request import urlopen
else:
  from urllib2 import urlopen

## at UWPKFCVM/parkfield_vptable.txt and parkfield_vstable.txt

model = "UWPKFCVM"

dimension_x = 0
dimension_y = 0 
dimension_z = 0

def usage():
    print("\n./make_data_files.py\n\n")
    sys.exit(0)

def download_urlfile(url,fname):
  print("\ndata file:",url,"\n")
  try:
    response = urlopen(url)
    CHUNK = 16 * 1024
    with open(fname, 'wb') as f:
      while True:
        chunk = response.read(CHUNK)
        if not chunk:
          break
        f.write(chunk)
  except:
    e = sys.exc_info()[0]
    print("Exception retrieving and saving model datafiles:",e)
    raise
  return True

def download_urlfile2(url, fname):
    print("\ndata file:", url, "\n")
    try:
        context = ssl.create_default_context(cafile=certifi.where())
        response = urlopen(url, context=context)
        CHUNK = 16 * 1024
        with open(fname, 'wb') as f:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                  break
                f.write(chunk)
    except Exception as e:
        print("Exception retrieving and saving model datafiles:", e)
        raise
    return True


def main():

    # Set our variable defaults.
    path = ""
    mdir = ""

    try:
        fp = open('./config','r')
    except:
        print("ERROR: failed to open config file")
        sys.exit(1)

    ## look for model_data_path and other varaibles
    lines = fp.readlines()
    for line in lines :
        if line[0] == '#' :
          continue
        parts = line.split('=')
        if len(parts) < 2 :
          continue;
        variable=parts[0].strip()
        val=parts[1].strip()

        if (variable == 'model_data_path') :
            path = val + '/' + model
            continue
        if (variable == 'model_dir') :
            mdir = "./"+val
            continue
        if (variable == 'nx') :
            dimension_x = int(val)
            continue
        if (variable == 'ny') :
            dimension_y = int(val)
            continue
        if (variable == 'nz') :
            dimension_z = int(val)
            continue
        continue
    if path == "" :
        print("ERROR: failed to find variables from config file")
        sys.exit(1)

    fp.close()

#    print("\nDownloading model file\n")
#
#    download_urlfile2(url,fname)

    subprocess.check_call(["mkdir", "-p", mdir])

    # Now we need to go through the data files and put them in the correct
    # format. More specifically, we need a vp.dat

    fvp = open("./uwpkfcvm_vp.txt", "r")
    fvs = open("./uwpkfcvm_vs.txt", "r")
    fvp_out = open("uwpkfcvm/vp.dat", "wb")
    fvs_out = open("uwpkfcvm/vs.dat", "wb")

    vp_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))
    vs_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))

    data_total_cnt=0
    x_pos=0
    y_pos=0
    z_pos=0
    for vpline, vsline in zip(fvp, fvs):
        in_vp = vpline.split()
        in_vs = vsline.split()

        skip_x = (in_vp[0])
        skip_y = (in_vp[1])
        val_z = float(in_vp[2])
        val_vp = float(in_vp[3])

        skip_x = (in_vs[0])
        skip_y = (in_vs[1])
        val_z = float(in_vs[2])
        val_vs = float(in_vs[3])

        val_vp = val_vp * 1000.0
        val_vs = val_vs * 1000.0
        val_z = val_z * 1000.0

        loc = z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
        vp_arr[loc]=val_vp
        vs_arr[loc]=val_vs

        if (data_total_cnt < 5) :
           print("%d : %lf %lf %lf\n" % (loc, val_vp, val_vs, val_z))
        data_total_cnt=data_total_cnt+1

        x_pos = x_pos + 1
        if(x_pos == dimension_x) :
          x_pos = 0;
          y_pos = y_pos+1
          if(y_pos == dimension_y) :
            y_pos=0;
            z_pos = z_pos+1
            if(z_pos == dimension_z) :
              print ("All DONE")


    vp_arr.tofile(fvp_out)
    vs_arr.tofile(fvs_out)

    fvp.close()
    fvs.close()
    fvp_out.close()
    fvs_out.close()
    print("Done! total data cnt = %d"%data_total_cnt)


if __name__ == "__main__":
    main()

