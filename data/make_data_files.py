#!/usr/bin/env python

##
#  Builds the data files in the expected format from 
#
# from >> XX??  YY?? depth(KM) Latitude Longitude  Vp/Vs (km/s)
#         -65.0 -70.0   0.0    35.10035 -120.55757 3.17
#
# with varying depth
#
## Origin and Rotate Angle(anticlockwise):35.960 -120.505 -40.00
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

    fvp = open("./parkfield_vptable.txt", "r")
    fvs = open("./parkfield_vstable.txt", "r")
    f_out = open("uwpkfcvm/parkfield.txt", "w")

    data_total_cnt=0
    vp_nan_cnt=0
    vs_nan_cnt=0

    for vpline, vsline in zip(fvp, fvs):
        arr_vp = vpline.split()
        arr_vs = vsline.split()

        skip_x = float(arr_vp[0])
        skip_y = float(arr_vp[1])
        in_z = float(arr_vp[2])
        in_lat = float(arr_vp[3])
        in_lon = float(arr_vp[4])

        tmp = arr_vp[5]
        if( tmp != "NaN" ) :
           in_vp = float(tmp)
           in_vp = in_vp * 1000.0;
        else:
           vp_nan_cnt = vp_nan_cnt + 1

        tmp = arr_vs[5]
        if( tmp != "NaN" ) :
           in_vs = float(tmp)
           in_vs = in_vs * 1000.0;
        else:
           vs_nan_cnt = vs_nan_cnt + 1

        data_total_cnt = data_total_cnt + 1

#-123.267 39.3481 0 763.174 2087.571 1957.221
        f_out.write("%lf %lf %lf %lf %lf\n" % (in_lon,in_lat,in_z,in_vs,in_vp))

    fvp.close()
    fvs.close()
    f_out.close()
    print("Done! total data cnt = %d"%data_total_cnt)


if __name__ == "__main__":
    main()

