#!/bin/bash

#Minimal runscript for atmospheric dynamical cores

#=====================================================================================================
# Parse command line options
#=====================================================================================================

Help()
{ echo "Summary of command line options:"
  echo "--------------------------------"
  echo "-s --skip-mkmf    Skip building of Makefile and go straight to "
  echo "                  compilation of model"
    }
TEMP=$(getopt -o s,h --long skip-mkmf,help -- "$@")

skip_mkmf=false

eval set -- "$TEMP"
while true ; do
    case "$1" in
	-s|--skip-mkmf       ) skip_mkmf=true ; shift 1;;
	-h|--help            ) Help ; exit 1 ; shift 1 ;;
	--                   ) shift ; break ;;
	*                    ) echo "Error parsing"; exit 1 ;;
    esac
done

#=====================================================================================================
# Source config -- ensure you have user_config file in this directory!
#=====================================================================================================
set -x 

if [[ ! -f user_config ]] ; then
    echo "No file named 'user_config' found. Create one with the same format as user_config_example"
    exit 1
fi
source user_config

#====================================================================================================
# Set important paths
#====================================================================================================

if [[ $openmp = true ]] ; then
    template=$fv3_home/tools/mkmf_templates/mkmf.template.${platform}_openmp   # path to mkmf template for platformx
else
    template=$fv3_home/tools/mkmf_templates/mkmf.template.$platform
fi
mkmf=$fv3_home/src/mkmf/bin/mkmf                                # path to executable mkmf
sourcedir=$fv3_home/src                                         # path to directory containing model source code
pp_path=$fv3_home/postprocessing/bin
run_script=$PWD/run_${run_name}              # path/name of this run script (for resubmit)
exp_home=$PWD                           # directory containing run/$run_script and input/
exp_name=${exp_home##*/}                       # name of experiment (i.e., name of this model build)
PATH="${PATH}:${pp_path}"

# Test that template actually exists
if [[ ! -f $template ]] ; then
    echo "mkmf template: $template does not exist!"
    echo "Make sure it can be found in $fv3_home/tools/mkmf_templates"
    exit 1
fi
#====================================================================================================
# Load required environment modules (edit for different systems)
#====================================================================================================
set +x
module load intel-compilers/2022
module load openmpi/4.1.5-intel
module load hdf5/1.14.0-intel-parallel
module load netcdf/netcdf-c-4.9.2-parallel
module load netcdf/netcdf-fortran-4.6.1
set -x

#====================================================================================================
# Setup working directory
#====================================================================================================
echo Working directory is $PWD
rm -rf workdir
# set initial conditions and move to executable directory
mkdir restart-files

if [[ $init_cond != "" ]] ; then 
    cp $init_cond/* $exp_home/restart-files
fi

execdir=$PWD/exec.$platform       # where code is compiled and executable is created
workdir=$PWD/workdir              # where model is run and model output is produced

pathnames=$PWD/path_names           # path to file containing list of source paths
namelist=$PWD/input.nml            # path to namelist file
diagtable=$PWD/diag_table           # path to diagnositics table
fieldtable=$PWD/field_table         # path to field table (specifies tracers)

if [[ ! -d $execdir ]] ; then mkdir $execdir ; fi
# if [[ -e $workdir ]] ; then
#   echo "ERROR: Existing workdir may contaminate run.  Move or remove $workdir and try again."
#   exit 1
# fi
mkdir $workdir $workdir/INPUT $workdir/RESTART
cp -r $cwd/restart-files/* $workdir/INPUT/

/bin/cp $template $workdir/tmp_template
echo "fv3_home = $fv3_home" >> $workdir/tmp_template

#====================================================================================================
# Find list of model source code and append to tmp_pathnames, to be read by mkmf
#====================================================================================================
# Prepend fortran files in srcmods directory to pathnames.
# Use 'find' to make list of srcmod/*.f90 files. mkmf uses only the first instance of any file name.

if [[ $skip_mkmf = false ]] ; then
    find $exp_home/srcmods/ -maxdepth 1 -iname "*.f90" -o -iname "*.F90" -o -iname "*.inc" -o -iname "*.c" -o -iname "*.h" > $workdir/tmp_pathnames
    echo "Using the following sourcecode modifications:"
    cat $workdir/tmp_pathnames

    set +x
    excludes=$fv3_home/tools/src_filter/exc_dirs # list of dirs to exclude from src search
    includes=$fv3_home/tools/src_filter/src_dirs # list of dirs to include from src search

    readarray includes < $includes
    readarray excludes < $excludes
    includes=$(eval echo ${includes[@]})
    excludes=($(eval echo ${excludes[@]}))

    exc_cmd=("( -path ${excludes[0]}")
    for exc in ${excludes[@]:1} ; do
	exc_cmd+=("-o -path $exc")
    done
    exc_cmd+=(")")

    find_args=($includes ${exc_cmd[@]} -prune -o \( -iname "*.f90" -o -iname "*.F90" -o -iname "*.inc" -o -iname "*.c" -o -iname "*.h" \) -print)
    set -x
    find "${find_args[@]}" >> $workdir/tmp_pathnames
fi

cd $execdir
#====================================================================================================
# Run mkmf to create makefile
#====================================================================================================

# compiler definitions for FMS (taken from autoconf build)

if [[ $skip_mkmf = false ]] ; then
    CDEFS_FMS=("-DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1\
       -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1\
       -DHAVE_DLFCN_H=1 -DLT_OBJDIR=\".libs/\" -DHAVE_MPI_H=1 -DHAVE_NETCDF_H=1 -DHAVE_GETTID=1\
       -DHAVE_SCHED_GETAFFINITY=1 -DHAVE_MOD_MPI=1 -DHAVE_MPIF_H=1 -DHAVE_MOD_NETCDF=1 \
       -DHAVE_CRAY_POINTER=1 -DHAVE_QUAD_PRECISION=1 -DHAVE_INTERNAL_NML=1 -Duse_netCDF=1 -Duse_libMPI=1")

    # compiler definitions for FV3
    CDEFS_FV3=("-DHI48 -DSPMD")

    CDEFS=("${CDEFS_FMS[@]} ${CDEFS_FV3[@]}")
    
    $mkmf -p fms.x -t $template -c "${CDEFS[@]}" -a $sourcedir\
	  $workdir/tmp_pathnames
fi

#====================================================================================================
# Compile the model
#====================================================================================================
make -j 4 -f Makefile
if [[ $? != 0 ]] ; then
  echo Compilation failed
  exit 1
fi

#====================================================================================================
# Run the model
#====================================================================================================
cd $workdir
cp $namelist input.nml
cp $diagtable diag_table
cp $fieldtable field_table
cp $execdir/fms.x fms.x
#--------------------------------------------------------------------------------------------------------
# run the model with mpirun
mpirun -np $npes --bind-to core --mca btl openib,self,vader --mca btl_openib_allow_ib 1 fms.x

if [[ $? != 0 ]] ; then echo "Error in model run" ; exit ; fi
#--------------------------------------------------------------------------------------------------------

#====================================================================================================
# Combine the netcdf file from each processor to one per grid tile (6 in total for cubedsphere)
#====================================================================================================
  if [[ $npes > 1 ]] ; then
     for ncfile in `/bin/ls *.nc.0000` ; do
  	mppnccombine -r -n4 ${ncfile%.*}
     done
  fi

# --------------------------------------------------------------------------------------------------------
if [[ $? != 0 ]] ; then  exit ; fi
#--------------------------------------------------------------------------------------------------------

#====================================================================================================
# Regrid data from each tile onto the latitude-longitude grid
#====================================================================================================
# Interpolate data to lat-lon grid
diagFiles=*.tile1.nc
latlonfiles=
mosaic=$fv3_home/postprocessing/mosaics/C$res/mosaic_C$res.nc


for File in $diagFiles ; do
  variables=`/usr/bin/ncdump -h $File | grep 'grid_yt, grid_xt' | awk '{print $2}' | cut -d\( -f1`
  variables=`echo $variables |sed 's/ /,/g'`
  basename=${File%.*.*}
  
  mpirun -np $npes fregrid_parallel --input_mosaic $mosaic --input_file $basename --interp_method conserve_order2\
	 --remap_file fregrid_remap_file --nlon $((res*4)) --nlat $((res*2)) --scalar_field $variables
   latlonfiles="$latlonfiles $basename.nc"
done
if [[ $? != 0 ]] ; then echo "Error in regrid" ; exit ; fi
echo 'Fields interpolated to lat-lon grid exist in these files:'

for File in $latlonfiles ; do
  ls -l $PWD/$File
done

# Ensure the new tiles have the hybrid sigma data for vertical regridding
ncks -A -v pk,bk atmos_static.tile1.nc atmos_daily.nc
ncks -A -v pk,bk atmos_static.tile1.nc atmos_average.nc

#=========================================================================================
# Vertically interpolate from hybrid sigma coord to constant pressure levels
#=========================================================================================
interp_files="$workdir/atmos_daily.nc $workdir/atmos_average.nc"


for File in $interp_files ; do
    pfull=$(ncdump -v pfull $File | sed -ne '/ pfull =/,$ p' | cut -f 2 -d '=' | cut -f 1 -d ';' | sed '$d' | sed 's/,/\ /g'| tr '\n' ' ')
    pfull=($pfull)

    # pfull_new=()

    # for p in "${pfull[@]}" ; do
    # 	test=$(printf %.0f $(echo "$p*100" | bc -l))

    # 	if [[ "$test" != 0 && ! "${pfull_new[*]}" =~ "$test" ]] ; then
    # 	    pfull_new=("${pfull_new[@]}" $test)
    # 	fi
    # done
    
    # pfull=("${pfull_new[@]}")
    (cd $pp_path ; plevel.bash -0 -a -p "${pfull[*]}" -i $File -o "${File%.*}_interp.nc")
done

# for File in $interp_files ; do
#     pfull=$(ncdump -v pfull $File | sed -ne '/ pfull =/,$ p' | cut -f 2 -d '=' | cut -f 1 -d ';' | sed '$d' | sed 's/,/\ /g'| tr '\n' ' ')
    
#     pfull_new=
#     set +x 
#     for p in $pfull ; do
# 	pfull_new="$pfull_new $(printf %.0f $(echo "$p*100"| bc -l) )"
#     done
#     set -x
#     pfull=$(echo $pfull_new | xargs)
#     (cd $pp_path ; plevel.bash -0 -a -p ''"$pfull"'' -i $File -o "${File%.*}_interp.nc")
# done
if [[ $? != 0 ]] ; then echo "Error in vertical interpolation" ; exit ; fi

#=========================================================================================
# Move output
#=========================================================================================
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

#=========================================================================================
# Plotting
#=========================================================================================
echo "Plotting"
plot_dir="$fv3_home/tools/plotting"

python $plot_dir/plot_temp.py $output_dir --plevs $plot_plevels

