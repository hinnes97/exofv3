#!/bin/bash

# Modules required for the postprocessing:
# - intel compilers
# - hdf5
# - netcdf
# - openmpi
# Ensures the correct environment variables are set for compilation
# Will be system dependent, change as necessary. These are set up for
# use on Oxford servers
#===================================================================
module load intel-compilers/2022
module load openmpi/4.1.5-intel
module load hdf5/1.14.0-intel-parallel
module load netcdf/netcdf-c-4.9.2-parallel
module load netcdf/netcdf-fortran-4.6.1
#===================================================================

nctools_dir="../src/FRE-NCtools"
install_dir="${PWD}/../postprocessing"
cd $nctools_dir

# Run autoreconf
autoreconf -i .
if [[ $? != 0 ]] ; then echo "Autoreconf failed, make sure this is available on your system" ; fi

# Make build directory
if [[ ! -d "build" ]] ; then
    mkdir build
fi
cd build

# Run configure
../configure --with-mpi=yes --with-netcdf=yes --with-netcdf-fortran=yes --prefix=$install_dir
if [[ $? != 0 ]] ; then echo "Configure failed" ; fi

# Run make
make
if [[ $? != 0 ]] ; then echo "Error in make" ; fi
make install
if [[ $? != 0 ]] ; then echo "Error in make install" ; fi
