#!/bin/bash

#Minimal runscript for atmospheric dynamical cores
set -x
#--------------------------------------------------------------------------------------------------------
# define variables
# Check fms_home has been passed as command line argument
#=====================================================================================================
# User edit this section only
#=====================================================================================================
fms_home=$net2/new_fms
platform=oxford_ubuntu_1804                                      # A unique identifier for your platform
output_root=$net2/output_newfms
init_cond=""                                                    # Directory of restart file
rst_time=
run_name=default_nonhydro
plot_plevels="50.0 500.0"                                       # Pressure levels to plot at, in mbar
#====================================================================================================

npes=24
template=$fms_home/mkmf/templates/mkmf.template.$platform   # path to mkmf template for platformx
mkmf=$fms_home/mkmf/bin/mkmf                      # path to executable mkmf
sourcedir=$fms_home/src                           # path to directory containing model source code
pp_path=$fms_home/postprocessing/bin

#MDH
#source /etc/csh.cshrc
# source ~/.cshrc
set +x
module load intel-compilers/2022
module load openmpi/4.1.5-intel
module load hdf5/1.14.0-intel-parallel
module load netcdf/netcdf-c-4.9.2-parallel
module load netcdf/netcdf-fortran-4.6.1
set -x
echo Working directory is $PWD

run_script=$PWD/run_${run_name}              # path/name of this run script (for resubmit)
exp_home=$PWD                           # directory containing run/$run_script and input/
exp_name=${exp_home##*/}                       # name of experiment (i.e., name of this model build)

${PATH}="${PATH}:${pp_path}"
###############
fms_run=$PWD
rm -rf workdir
# set initial conditions and move to executable directory
mkdir restart-files

if [[ $init_cond != "" ]] ; then 
    cp $init_cond/* $fms_run/restart-files
fi

##############

#mppnccombine=$fms_home/FRE-NCtools/postprocessing/mppnccombine/mppnccombine

#--------------------------------------------------------------------------------------------------------
execdir=$PWD/exec.$platform       # where code is compiled and executable is created
workdir=$PWD/workdir              # where model is run and model output is produced

pathnames=$PWD/path_names           # path to file containing list of source paths
namelist=$PWD/input.nml            # path to namelist file
diagtable=$PWD/diag_table           # path to diagnositics table
fieldtable=$PWD/field_table         # path to field table (specifies tracers)

#--------------------------------------------------------------------------------------------------------
# setup directory structure
if [[ ! -d $execdir ]] ; then mkdir $execdir ; fi
if [[ -e $workdir ]] ; then
  echo "ERROR: Existing workdir may contaminate run.  Move or remove $workdir and try again."
  exit 1
fi
mkdir $workdir $workdir/INPUT $workdir/RESTART

#copy restart files
cp -r $cwd/restart-files/* $workdir/INPUT/

#--------------------------------------------------------------------------------------------------------
# compile the model code and create executable

# MDH
# append fms_home (containing netcdf libraries and include files) to template
/bin/cp $template $workdir/tmp_template
echo "fms_home = $fms_home" >> $workdir/tmp_template

# Prepend fortran files in srcmods directory to pathnames.
# Use 'find' to make list of srcmod/*.f90 files. mkmf uses only the first instance of any file name.
cd $sourcedir
find $exp_home/srcmods/ -maxdepth 1 -iname "*.f90" -o -iname "*.F90" -o -iname "*.inc" -o -iname "*.c" -o -iname "*.h" > $workdir/tmp_pathnames
echo "Using the following sourcecode modifications:"
cat $workdir/tmp_pathnames

find $sourcedir/atmos_exo -iname "*.f90" -o -iname "*.F90" -o -iname "*.c" -o -iname "*.h" -o -iname "*.inc" >> $workdir/tmp_pathnames
find $sourcedir -path $sourcedir/atmos_exo -prune -o \( -iname "*.f90" -o -iname "*.F90" -o -iname "*.c" -o -iname "*.h" -o -iname "*.inc" \) -print >> $workdir/tmp_pathnames

cd $execdir

$mkmf -p fms.x -t $template -c "-DHI_48 -Duse_libMPI -Duse_netCDF -DSPMD" -a $sourcedir $workdir/tmp_pathnames $workdir/tmp_pathnames_soc_int $workdir/tmp_pathnames_soc_src $workdir/tmp_pathnames_soc_mod $sourcedir/shared/mpp/include $sourcedir/shared/include
make -j 4 -f Makefile

if [[ $? != 0 ]] ; then
  echo Compilation failed
  exit 1
fi

cd $workdir
cp $namelist input.nml
cp $diagtable diag_table
cp $fieldtable field_table
cp $execdir/fms.x fms.x
#--------------------------------------------------------------------------------------------------------
# run the model with mpirun
mpirun -np $npes fms.x
if [[ $? != 0 ]] ; then echo "Error in model run" ; exit ; fi
#--------------------------------------------------------------------------------------------------------
#combine netcdf files
 if [[ $npes > 1 ]] ; then
    for ncfile in `/bin/ls *.nc.0000` ; do
 	mppnccombine -r -n4 ${ncfile%.*}
    done
 fi
#
#   --------------------------------------------------------------------------------------------------------
if [[ $? != 0 ]] ; then  exit ; fi
#--------------------------------------------------------------------------------------------------------
# Interpolate data to lat-lon grid
diagFiles=*.tile1.nc
latlonfiles=
mosaic_dir=$fms_home/FRE-NCtools/tools/make_solo_mosaic/C24

#cp $mosaic_dir/horizontal_grid.tile?.nc ./
#cp $mosaic_dir/mosaic_n48.nc ./

#fregrid=$fms_home/FRE-NCtools/tools/fregrid/fregrid_parallel

CXX=24
CXX2=48
CXX4=96

for File in $diagFiles ; do
  variables=`/usr/bin/ncdump -h $File | grep 'grid_yt, grid_xt' | awk '{print $2}' | cut -d\( -f1`
  variables=`echo $variables |sed 's/ /,/g'`
  basename=${File%.*.*}
  
   mpirun -np $npes fregrid_parallel --input_mosaic $mosaic_dir/mosaic_n48.nc --input_file $basename --interp_method conserve_order2 --remap_file fregrid_remap_file --nlon $CXX4 --nlat $CXX2 --scalar_field $variables 
   latlonfiles="$latlonfiles $basename.nc"
done
if [[ $? != 0 ]] ; then echo "Error in regrid" ; exit ; fi
echo 'Fields interpolated to lat-lon grid exist in these files:'

for File in $latlonfiles ; do
  ls -l $PWD/$File
done

ncks -A -v pk,bk atmos_static.tile1.nc atmos_daily.nc
ncks -A -v pk,bk atmos_static.tile1.nc atmos_average.nc

interp_files="$workdir/atmos_daily.nc $workdir/atmos_average.nc"
#interp_dir="$fms_home/FRE-NCtools/postprocessing/plevel"
#interp_script="plevel.sh"
for File in $interp_files ; do
    pfull=$(ncdump -v pfull $File | sed -ne '/ pfull =/,$ p' | cut -f 2 -d '=' | cut -f 1 -d ';' | sed '$d' | sed 's/,/\ /g'| tr '\n' ' ')
    pfull_new=
    set +x 
    for p in $pfull ; do
	pfull_new="$pfull_new $(printf %.0f $(echo "$p*100"| bc -l) )"
    done
    set -x
    pfull=$(echo $pfull_new | xargs)
    (cd $pp_path ; plevel.sh -0 -a -p ''"$pfull"'' -i $File -o "${File%.*}_interp.nc")
done
if [[ $? != 0 ]] ; then echo "Error in vertical interpolation" ; exit ; fi

# Move output
# Get runtime
runtime=$(grep 'days' input.nml | awk '{print $3}')
output_dir="$output_root/$run_name/$((runtime + rst_time))"

if [[ -d $output_dir ]] ; then  rm -rf $output_dir ; fi
mkdir -p $output_dir
rm -f atmos*tile*.nc
mv atmos_* $output_dir
cp input.nml $output_dir
cp ../run.bash $output_dir
mv RESTART $output_dir
mv logfile* $output_dir

cp -rf srcmods $output_dir

echo "Moved data" 

echo "Plotting"
conda activate analyse
plot_dir="$fms_home/plotting"

python $plot_dir/plot_temp.py $output_dir --plevs $plot_plevels

