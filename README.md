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
- The OpenMPI library, version 1.6 or newer should be installed, including its headers. Earlier versions have a bug limiting the number of concurrent processes to 128.
- Choose a directory in your host computer, e.g. `prato`

	`$> mkdir prato`
	`$> cd prato`

- Download the performance metrics library from https://github.com/lichinka/performance_metrics (you need the PAPI library from http://icl.cs.utk.edu/papi/software/index.html):

	`prato> git clone https://github.com/lichinka/performance_metrics performance`
	`prato> cd performance`
	`prato/performance> make clean && make`
	`prato/performance> cd ..`

- Download the OpenCL common library from https://github.com/lichinka/ocl_common You need OpenCL working on your computer before doing this. If you don't plan to use a GPU, you may download the AMD OpenCL SDK for your architecture from http://developer.amd.com/tools/heterogeneous-computing/amd-accelerated-parallel-processing-app-sdk/downloads/ and install it accordingly. After that:

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

- If anything goes wrong, make sure you have all the required libraries and that the variables in `Makefile.inc` point to the correct locations. On the other hand, if using MPICH2, you will have to adjust the `mpicc` command, for example:

	`$(shell mpicc --showme:compile)` => `-I/usr/include/mpich2`
	`$(shell mpicc --showme:link)` => `-L/usr/lib/mpich2`
	


Using PRATO
===========
- Start the GRASS environment, selecting or creating a `grassdata` directory where all GIS files will be saved.
- Adjust the dynamic library path so that it will find all the required libraries:

	`grass> export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:prato/src`

where `prato/src` is the directory containing the module binary and required libraries. 
- You may check if all the needed libraries are found at runtime with `ldd`. For example:

	`grass> ldd r.coverage`
	linux-vdso.so.1 =>  (0x00007fff275ff000)
	libgrass_gis.so => /usr/local/src/grass-6.4.2/dist.x86_64-unknown-linux-gnu/lib/libgrass_gis.so (0x00007ffec1a32000)
	libgrass_datetime.so => /usr/local/src/grass-6.4.2/dist.x86_64-unknown-linux-gnu/lib/libgrass_datetime.so (0x00007ffec1827000)
	libgrass_dbmibase.so => /usr/local/src/grass-6.4.2/dist.x86_64-unknown-linux-gnu/lib/libgrass_dbmibase.so (0x00007ffec1614000)
	libgrass_dbmiclient.so => /usr/local/src/grass-6.4.2/dist.x86_64-unknown-linux-gnu/lib/libgrass_dbmiclient.so (0x00007ffec1408000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007ffec115b000)
	libmpich.so.3 => /usr/lib/libmpich.so.3 (0x00007ffec0d7f000)
	libmpl.so.1 => /usr/lib/libmpl.so.1 (0x00007ffec0b7a000)
	libperformance.so => /home/grassuser/prato/src/libperformance.so (0x00007ffec0572000)
	libworker.so => /home/grassuser/prato/src/libworker.so (0x00007ffec0366000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007ffebffc3000)
	libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007ffebfdaa000)
	libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007ffebfba6000)
	librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007ffebf99d000)
	libcr.so.0 => /usr/lib/libcr.so.0 (0x00007ffebf793000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007ffebf576000)
	/lib64/ld-linux-x86-64.so.2 (0x00007ffec1c9f000)
	libpapi.so => /home/grassuser/papi/lib/libpapi.so (0x00007ffebf327000)
	liboclcommon.so => /home/grassuser/prato/src/liboclcommon.so (0x00007ffebf122000)
	libOpenCL.so.1 => /home/grassuser/AMD_SDK/lib/x86_64/libOpenCL.so.1 (0x00007ffebef1c000)
	libpfm.so.4 => /home/grassuser/papi/lib/libpfm.so.4 (0x00007ffebec2e000)

- Open the `parameters.ini` file and set the values according to your installation. Start by changing the DEM and Clutter map names. 
- It is important to keep all GIS layers (e.g. DEM, clutter, ...) in the same resolution! Check them by running:

	`grass> g.region -p | grep res`
	`nsres:      100`
	`ewres:      100`
