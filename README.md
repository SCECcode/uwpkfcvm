# UWPKFCVM (Parkfield)

<a href="https://github.com/sceccode/uwpkfcvm.git"><img src="https://github.com/sceccode/uwpkfcvm/wiki/images/uwpkfcvm_logo.png"></a>

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/uwpkfcvm)
[![uwpkfcvm-ucvm-ci Actions Status](https://github.com/SCECcode/uwpkfcvm/workflows/uwpkfcvm-ucvm-ci/badge.svg)](https://github.com/SCECcode/uwpkfcvm/actions)


 The UW Parkfield velocity model describes 3D seismic P- and S-wave velocities from a tomographic 
inversion for crustal structure in the greater Parkfield region (Zeng et al., 2016). The 
inversion used the code tomoDD by Zhang and Thurber (BSSA, 2003), which inverts body-wave 
arrival times for 3D P- and S-wave velocity structure simultaneously with determining 
earthquake locations. A novel aspect of this model is that it used P- and S-wave arrivals 
from low-frequency earthquakes (LFEs) to improve resolution at greater depths than previous 
body-wave tomography studies of the Parkfield region. The coordinate system origin is 
(35° 57.60' N, -120° 30.28' W) with the Y-axis rotation to an azimuth of 139.2°.

Zeng, Xiangfang, C. H. Thurber, D. R. Shelly, R. M. Harrington, E. S. Cochran, N. L. Bennington, 
D. Peterson, B. Guo and K. McClement, '3-D P- and S-wave velocity structure and low-frequency 
earthquake locations in the Parkfield, California region'.  Geophys. J. Int. (2016) 206, 1574–1585.
doi: 10.1093/gji/ggw217

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.7 or higher.

<pre>
aclocal -I m4
autoconf
automake --add-missing --force-missing

./configure --prefix=$UCVM_INSTALL_PATH/model/uwpkfcvm --with-proj-libdir=$UCVM_INSTALL_PATH/lib/proj/lib --with-proj-incdir=$UCVM_INSTALL_PATH/lib/proj/include --with-tiff-libdir=$UCVM_INSTALL_PATH/lib/tiff/lib --with-tiff-incdir=$UCVM_INSTALL_PATH/lib/tiff/include --with-sqlite-libdir=$UCVM_INSTALL_PATH/lib/sqlite/lib --with-sqlite-incdir=$UCVM_INSTALL_PATH/lib/sqlite/include
</pre>

## Contact the authors

If you would like to contact the authors regarding this software,
please e-mail software@scec.org. Note this e-mail address should
be used for questions regarding the software itself (e.g. how
do I link the library properly?). Questions regarding the model's
science (e.g. on what paper is the UWPKFCVM based?) should be directed
to the model's authors, located in the AUTHORS file.

## NOTE

A right rectangle with ??? degree rotation

density are calculated, ??? ....


