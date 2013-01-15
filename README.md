A high-performance parallel radio coverage prediction tool for GRASS GIS.-


Installing GRASS
================
- Download GRASS source code from http://grass.fbk.eu/download/software.php#g64x
- Unpack the .tar.gz file in /usr/local/src (as root).
- For freshly installed machines, you will need PROJ4, GDAL, FFTW, Python 2.x and some GUI toolkit (WxWidgets, Tcl/Tk, ...) installed. The method on how to provide these dependencies differs among distributions. To indicate a different Python interpreter, be sure that a symlink exists to the correct one, e.g.:
    
	`$> ls -lah /usr/bin/python`
	`lrwxrwxrwx 1 root root 16 Feb 11 11:00 /usr/bin/python -> /usr/bin/python2`

or change the PYTHON setting in include/Make/Platform.make and Platform.make.in, e.g.:
    
	`PYTHON=python2`

- Compile and install GRASS from the source code (as root). Binary installations (i.e. from a package manager) are not good for GRASS module development. Even the `grass-6.x.x-dev` packages under Debian or Ubuntu-like distributions have bugs, so do not use them. 

- Example to install GRASS from source (adjust options for your machine/needs):

	`$> ./configure --with-sqlite \
			--enable-64bit \
			--with-tcltk \
			--with-python=/usr/bin/python2.7-config \
			--prefix=/usr/local`

- After `configure` has finished correctly, issue `make` and the `make install`. The resulting installation should now be located at `/usr/local/grass-6.4xxx`, including the binaries (e.g. `grass` and/or `grass64`) to start the program, located at `/usr/local/bin`.
- If you are planning to do GRASS module development, change both directories ownerships (i.e. the source code tree and the installation one) to avoid module compilation as root. For example, if your everyday user is a member of the `users` group, do:

	`$> sudo chown --recursive root:users /usr/local/grass-6.4xxx`
	`$> sudo chown --recursive root:users /usr/local/src/grass-6.4xxx`
	`$> sudo chmod --recursive g+w /usr/local/grass-6.4xxx`
	`$> sudo chmod --recursive g+w /usr/local/src/grass-6.4xxx`


Installing PRATO
================

Requirements
------------
- The OpenMPI library, version 1.6 or newer should be installed, including its headers.
- Choose a directory in your host computer, e.g. `prato`

	`$> mkdir prato`
	`$> cd prato`

- Download the performance metrics library from https://github.com/lichinka/performance_metrics (you need the PAPI library):

	`prato> git clone https://github.com/lichinka/performance_metrics performance`
	`prato> cd performance`
	`prato/performance> make clean && make`
	`prato/performance> cd ..`

- Download the OpenCL common library from https://github.com/lichinka/ocl_common (you need OpenCL working on your computer before doing this):

	`prato> git clone https://github.com/lichinka/ocl_common ocl_common`
	`prato> cd ocl_common`
	`prato/ocl_common> make clean && make && ./cl_test`
	`prato/ocl_common> cd ..`

- Download the lastest version of PRATO from https://github.com/lichinka/prato 

	`prato> git clone https://github.com/lichinka/prato src`

- Create symlinks of the dependency libraries into the `src` directory. This will ease the dynamic-library path handling at runtime:

	`prato> cd src`
	`prato/src> ln -s ../performance/libperformance.so .`
	`prato/src> ln -s ../ocl_common/liboclcommon.so .`


Compilation
-----------
- The dir structure within `prato` should be:

	`prato`
	`  |`
	`  |--- performance`
	`  |--- ocl_common`
	`  |--- src`

where `performance` and `ocl_common` are the dependency libraries listed above.

- Compile the source code of PRATO to create a GRASS module:

	`prato/src> make clean && make`

- If anything goes wrong, make sure you have all the required libraries and that the variables in `Makefile.inc` point to the correct locations.


Using PRATO
===========
- Start the GRASS environment, selecting or creating a `grassdata` directory where all GIS files will be saved.
- It is important to keep all GIS layers (e.g. DEM, clutter, ...) in the same resolution! Check them by running:

	`grass> g.region -p | grep res`
	`nsres:      100`
	`ewres:      100`

