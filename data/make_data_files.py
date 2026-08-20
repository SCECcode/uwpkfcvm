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
import shutil
import subprocess
import os
import struct
import array

model = "UWPKFCVM"

dimension_x = 0
dimension_y = 0 
dimension_z = 0

def usage():
    print("\n./make_data_files.py\n\n")
    sys.exit(0)

def download_urlfile(url, fname):
    # Option 1A: aria2c tuned for slow/unstable connections
    if shutil.which("aria2c"):
        cmd = [
            "aria2c",
            "-x", "4",               # Limit to 4 connections (prevents network congestion)
            "-s", "4",               # Split into 4 parts
            "-c",                    # Always resume partial downloads
            "--max-tries=0",         # Infinite retries if Wi-Fi drops
            "--retry-wait=5",        # Wait 5 sec between retries
            "-o", fname,
            url
        ]
    # Option 1B: curl with resume fallback
    elif shutil.which("curl"):
        cmd = [
            "curl",
            "-L",                    # Follow redirects
            "-C", "-",               # Resume automatically
            "--retry", "999",        # Retry on failure
            "--retry-delay", "5",
            "-o", fname,
            url
        ]
    else:
        raise RuntimeError("Neither aria2c nor curl is installed.")

    process = subprocess.run(cmd, check=True)
# Check for success
    if process.returncode == 0 and os.path.exists(fname):
        print(f"\n[SUCCESS] Download completed! Proceeding with script...")
        return True
    else:
        raise RuntimeError(f"Download failed with exit code {process.returncode}")
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

## at UWPKFCVM/uwpkfcvm_vp.txt and uwpkfcvm_vs.txt
    print("\nDownloading model file\n")

    url=path+"/uwpkfcvm_vp.txt";
    download_urlfile(url,"./uwpkfcvm_vp.txt")
    url=path+"/uwpkfcvm_vs.txt";
    download_urlfile(url,"uwpkfcvm_vs.txt")

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

